"""
Checks for the geometry in pass_predictor.py. These do not need SGP4.

Run with:  python test_geometry.py
"""

import numpy as np

from pass_predictor import (compass, find_passes, gmst, julian_date,
                            look_angles, observer_ecef, parse_tle, teme_to_ecef)


def test_julian_date_at_j2000():
    # 2000-01-01 12:00:00 UTC is Julian date 2451545.0
    jd, fr = julian_date(946728000)
    assert abs(jd + fr - 2451545.0) < 1e-9


def test_gmst_at_j2000():
    # Reference value: 280.4606 degrees
    jd, fr = julian_date(946728000)
    assert abs(np.degrees(gmst(jd, fr)) - 280.4606) < 1e-3


def test_observer_on_equator():
    x, y, z = observer_ecef(0.0, 0.0, 0.0)
    assert abs(x - 6378.137) < 1e-6 and abs(y) < 1e-9 and abs(z) < 1e-9


def test_observer_at_pole():
    x, y, z = observer_ecef(90.0, 0.0, 0.0)
    assert abs(z - 6356.752) < 1e-3


def test_satellite_straight_overhead():
    lat, lon = 21.2703, 40.4158
    sat = observer_ecef(lat, lon, 400.0).reshape(1, 3)
    az, el, rng = look_angles(sat, lat, lon, 0.0)
    assert abs(el[0] - 90.0) < 1e-4
    assert abs(rng[0] - 400.0) < 1e-6


def test_satellite_due_north_on_horizon_side():
    sat = observer_ecef(10.0, 0.0, 400.0).reshape(1, 3)
    az, el, _ = look_angles(sat, 0.0, 0.0, 0.0)
    assert abs(az[0]) < 1e-6 or abs(az[0] - 360.0) < 1e-6
    assert 0 < el[0] < 90


def test_teme_to_ecef_quarter_turn():
    r = teme_to_ecef(np.array([[1.0, 0.0, 0.0]]), np.array([np.pi / 2]))
    assert np.allclose(r, [[0.0, -1.0, 0.0]], atol=1e-12)


def test_find_passes_ignores_cut_off_pass():
    t = np.arange(10.0)
    el = np.array([-5, 5, 15, 30, 15, 5, -5, 12, 20, 25], dtype=float)
    passes = find_passes(t, np.zeros(10), el, 10.0)
    assert len(passes) == 1
    assert passes[0]["peak_el"] == 30 and passes[0]["rise"] == 2 and passes[0]["set"] == 4


def test_compass():
    assert [compass(a) for a in (0, 44, 90, 180, 270, 350)] == ["N", "NE", "E", "S", "W", "N"]


def test_parse_tle_with_name():
    text = ("ISS (ZARYA)\n"
            "1 25544U 98067A   26148.13113954  .00011691  00000+0  21663-3 0  9996\n"
            "2 25544  51.6335  39.3887 0007375 106.1024 254.0777 15.49434162568649\n")
    name, line1, line2 = parse_tle(text)
    assert name == "ISS (ZARYA)" and line1.startswith("1 25544") and line2.startswith("2 25544")


if __name__ == "__main__":
    tests = [f for name, f in sorted(globals().items()) if name.startswith("test_")]
    for test in tests:
        test()
        print(f"ok  {test.__name__}")
    print(f"\n{len(tests)} tests passed")
