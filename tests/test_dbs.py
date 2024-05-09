import pickle

import numpy
import pytest

from cogent3 import load_table

from ensembl_lite._aligndb import AlignDb, AlignRecord
from ensembl_lite._homologydb import HomologyDb
from ensembl_lite._install import LoadHomologies, _load_one_align


@pytest.fixture(scope="function")
def db_align(DATA_DIR, tmp_path):
    records = _load_one_align()(DATA_DIR / "tiny.maf")
    outpath = tmp_path / "blah.sqlitedb"
    db = AlignDb(source=outpath)
    db.add_records(records)
    db.close()
    return AlignDb(source=outpath)


def test_db_align(db_align):
    orig = len(db_align)
    source = db_align.source
    db_align.close()
    got = AlignDb(source=source)
    assert len(got) == orig


def test_db_align_add_records(db_align):
    gap_spans = numpy.array([[2, 5], [7, 1]], dtype=int)
    orig = dict(
        source="blah",
        block_id="42",
        species="human",
        seqid="1",
        start=22,
        stop=42,
        strand="-",
    )
    records = [AlignRecord(gap_spans=gap_spans, **orig.copy())]
    db_align.add_records(records)
    got = list(
        db_align._execute_sql(
            f"SELECT * from {db_align.table_name} where block_id=?", (42,)
        )
    )[0]
    index = got["id"]
    got = {k: got[k] for k in orig if k != "id"}
    # gap_spans not stored in the db, but in a GapStore
    assert got == orig
    assert index == 1


@pytest.mark.parametrize("func", (str, repr))
def test_db_align_repr(db_align, func):
    func(db_align)


@pytest.fixture
def hom_dir(DATA_DIR, tmp_path):
    path = DATA_DIR / "small_protein_homologies.tsv.gz"
    table = load_table(path)
    outpath = tmp_path / "small_1.tsv.gz"
    table[:1].write(outpath)
    outpath = tmp_path / "small_2.tsv.gz"
    table[1:2].write(outpath)
    return tmp_path


def test_extract_homology_data(hom_dir):
    loader = LoadHomologies(
        {"gorilla_gorilla", "nomascus_leucogenys", "notamacropus_eugenii"}
    )
    got = loader(hom_dir.glob("*.tsv.gz"))
    assert len(got) == 2
    # loader dest cols matches the db schema
    assert set(loader.dest_col) == HomologyDb._homology_schema.keys()


def test_homology_db(hom_dir):
    loader = LoadHomologies(
        {"gorilla_gorilla", "nomascus_leucogenys", "notamacropus_eugenii"}
    )
    got = loader(hom_dir.glob("*.tsv.gz"))
    outpath = hom_dir / "species.sqlitedb"
    db = HomologyDb(source=outpath)
    db.add_records(records=got, col_order=loader.dest_col)
    assert len(db) == 2
    db.close()
    got = HomologyDb(source=outpath)
    assert len(got) == 2


def test_pickling_db(db_align):
    # should not fail
    pkl = pickle.dumps(db_align)
    upkl = pickle.loads(pkl)
    assert db_align.source == upkl.source
