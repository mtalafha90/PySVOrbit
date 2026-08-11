import backend
import pysvorbit


def test_public_import_facade_uses_backend_engine():
    assert pysvorbit.orb is backend.orb
    assert pysvorbit.readinp is backend.readinp
    assert pysvorbit.fitorb is backend.fitorb
    assert pysvorbit.eph is backend.eph


def test_public_version_matches_distribution_release():
    assert pysvorbit.__version__ == "1.0.1"
