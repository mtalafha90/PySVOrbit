"""Exercise independently referenced Besselian, JD and MJD conversions."""
import numpy as np
import pytest
import backend

# USNO Circular 163, Appendix: B1900 = JD 2415020.31352 and
# one Besselian year = 365.242198781 days.
JD_B1900 = 2415020.31352
DAYS_PER_BESSELIAN_YEAR = 365.242198781

@pytest.mark.parametrize('epoch', [1900.0, 1950.0, 2000.0, 2015.1234])
def test_besselian_to_mjd_and_jd(epoch):
    jd = JD_B1900 + DAYS_PER_BESSELIAN_YEAR * (epoch - 1900.0)
    for t0, expected in [(57000.0, jd - 2400000.5), (2457000.0, jd)]:
        observations = np.array([[epoch, 1.0, 0.1]], dtype=float)
        backend.correct(observations, t0)
        assert observations[0, 0] == pytest.approx(expected, rel=0, abs=1e-8)

@pytest.mark.parametrize('epoch', [1900.0, 1950.0, 2000.0, 2015.1234])
def test_mjd_and_jd_to_besselian(epoch):
    jd = JD_B1900 + DAYS_PER_BESSELIAN_YEAR * (epoch - 1900.0)
    for observation_epoch in (jd, jd - 2400000.5):
        observations = np.array([[observation_epoch, 1.0, 0.1]], dtype=float)
        backend.correct(observations, 2000.0)
        assert observations[0, 0] == pytest.approx(epoch, rel=0, abs=1e-8)

def test_mjd_and_jd_conversion_does_not_shift_half_day():
    observations = np.array([[57000.0, 1, 1], [2457000.5, 1, 1]], dtype=float)
    backend.correct(observations, 2457000.0)
    assert observations[:, 0] == pytest.approx([2457000.5, 2457000.5])
