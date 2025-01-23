# consumes homology data from Ensembl tsv files and writes a parquet file
# for representation by _homology.HomologyDb
import pathlib
import shutil
import typing

import duckdb
from cogent3.app import typing as c3_types
from cogent3.app.composable import LOADER, define_app

from ensembl_tui import _config as eti_config
from ensembl_tui import _homology as eti_homology
from ensembl_tui import _ingest_annotation as eti_annotation
from ensembl_tui import _util as eti_util

T = dict[str, tuple[eti_homology.homolog_group, ...]]


def grouped_related(
    data: typing.Iterable[tuple[str, str, str, str, str]],
) -> T:
    """determines related groups of genes

    Parameters
    ----------
    data
        list of [(rel_type, sp1, gid1, sp2, gid2), ...]

    Returns
    -------
    a data structure that can be json serialised

    Notes
    -----
    I assume that for a specific relationship type, a gene can only belong
    to one group.
    """
    # grouped is {<relationship type>: {gene id: homolog_group}. So gene's
    # that belong to the same group have the same value
    grouped = {}
    for rel_type, sp1, gene_id_1, sp2, gene_id_2 in data:
        relationship = grouped.get(rel_type, {})
        if gene_id_1 in relationship:
            val = relationship[gene_id_1]
        elif gene_id_2 in relationship:
            val = relationship[gene_id_2]
        else:
            val = eti_homology.homolog_group(relationship=rel_type)
        val.gene_ids |= {gene_id_1: sp1, gene_id_2: sp2}

        relationship[gene_id_1] = relationship[gene_id_2] = val
        grouped[rel_type] = relationship

    return {
        rel_type: tuple(set(groups.values())) for rel_type, groups in grouped.items()
    }


def merge_grouped(
    groups: typing.Sequence[eti_homology.homolog_group],
) -> tuple[eti_homology.homolog_group, ...]:
    """merges grouped homology data"""
    grouped = {}
    for group in groups:
        if present := group.gene_ids.keys() & grouped.keys():
            val = grouped[next(iter(present))]
            val.gene_ids |= group.gene_ids
        else:
            grouped |= {gene_id: group for gene_id in group.gene_ids}

    return tuple(set(grouped.values()))


@define_app(app_type=LOADER)
class load_homologies:  # noqa: N801
    """app to load homology groups from a single Ensembl tsv file"""

    def __init__(self, allowed_species: set[str]) -> None:
        # map the Ensembl columns to HomologyDb columns
        self._src_cols = (
            "homology_type",
            "species",
            "gene_stable_id",
            "homology_species",
            "homology_gene_stable_id",
        )
        self._create_sql = (
            f"CREATE TABLE my_table AS SELECT {','.join(self._src_cols)} FROM"
            " read_csv_auto('{}', delim='\t', header=True)"
        )
        allowed = ", ".join(f"{sp!r}" for sp in allowed_species)
        self._select_sql = (
            f"SELECT {','.join(self._src_cols)} FROM my_table"
            f" WHERE species IN ({allowed}) AND homology_species IN ({allowed})"
        )

    def main(self, path: c3_types.IdentifierType) -> c3_types.SerialisableType:
        conn = duckdb.connect(":memory:")
        sql = self._create_sql.format(path)
        conn.sql(sql)
        return grouped_related(conn.sql(self._select_sql).fetchall())


class HomologyAggregator:
    def __init__(self, conn: duckdb.DuckDBPyConnection) -> None:
        self._homology_id = 0
        self.conn = conn

    def add_records(
        self,
        *,
        records: typing.Sequence[eti_homology.homolog_group],
        relationship_type: str,
    ) -> None:
        """inserts homology data from records

        Parameters
        ----------
        records
            a sequence of homolog group instances, all with the same
            relationship type
        relationship_type
            the relationship type
        """
        if not relationship_type:
            msg = f"invalid {relationship_type=!r}"
            raise ValueError(msg)

        values = []
        for homology_id, group in enumerate(records, start=self._homology_id + 1):
            if group.relationship != relationship_type:
                msg = f"{group.relationship=} != {relationship_type=}"
                raise ValueError(msg)

            # get geneids and species for this group, storing the
            # geneid id for each record

            values.extend(
                [
                    (homology_id, gene_id, species, relationship_type)
                    for gene_id, species in group.gene_ids.items()
                ],
            )

        # record last homology_index
        self._homology_id = homology_id
        # create the homology table entries
        placeholder = ", ".join("?" * (len(eti_homology.HOMOLOGY_ATTR_COLS) - 1))
        cols = ", ".join(eti_homology.HOMOLOGY_ATTR_COLS[1:])
        sql = (
            f"INSERT OR IGNORE INTO homology_groups_attr({cols}) VALUES ({placeholder})"
        )
        self.conn.executemany(sql, parameters=values)


def make_homology_aggregator_db() -> HomologyAggregator:
    conn = duckdb.connect(":memory:")
    # this view condenses the homology member data so that we can
    # query for a relationship type and get the stableid and species
    # bear in mind that for a given relationship type, I'm assuming
    # that a stableid can only belong to one homology group (homology_id)
    conn.sql("CREATE SEQUENCE rowid_seq")
    schema = ", ".join(eti_homology.HOMOLOGY_ATTR_SCHEMA)
    sql = f"CREATE TABLE homology_groups_attr ({schema})"
    conn.sql(sql)
    return HomologyAggregator(conn)


def write_homology_views(agg: HomologyAggregator, outdir: pathlib.Path) -> None:
    # we write out the gene_species_attr and homology_groups_attr views
    # to the output directory
    outdir.mkdir(parents=True, exist_ok=True)
    _ = eti_annotation.export_parquet(
        con=agg.conn,
        table_name="homology_groups_attr",
        dest_dir=outdir,
    )


def local_install_homology(
    config: eti_config.Config,
    force_overwrite: bool = False,
    max_workers: int = 1,
    verbose: bool = False,
):
    from rich.progress import track

    if force_overwrite:
        shutil.rmtree(config.install_homologies, ignore_errors=True)

    config.install_homologies.mkdir(parents=True, exist_ok=True)

    tsv_paths = []
    for sp in config.db_names:
        path = config.staging_homologies / sp
        tsv_paths.extend(list(path.glob("*.tsv*")))

    if verbose:
        eti_util.print_colour(f"homologies {max_workers=}", "yellow")

    loader = load_homologies(allowed_species=set(config.db_names))

    tasks = eti_util.get_iterable_tasks(
        func=loader,
        series=tsv_paths,
        max_workers=max_workers,
    )
    agg = make_homology_aggregator_db()
    for grouped in track(tasks, total=len(tsv_paths), description="aggregating"):
        for rel_type, records in grouped.items():
            agg.add_records(records=records, relationship_type=rel_type)
