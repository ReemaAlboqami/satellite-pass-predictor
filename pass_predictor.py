"""
Satellite Pass Predictor

Predicts when a satellite passes over a ground location, using real
orbital data (TLEs) from CelesTrak and the SGP4 propagator.

Example:
    python pass_predictor.py                      # ISS over Taif, next 3 days
    python pass_predictor.py --catnr 20580        # Hubble Space Telescope
    python pass_predictor.py --lat 24.71 --lon 46.68 --alt 612   # Riyadh
"""

import argparse
import urllib.request
from datetime import datetime, timedelta, timezone
from pathlib import Path

import numpy as np

CELESTRAK_URL = "https://celestrak.org/NORAD/elements/gp.php?CATNR={catnr}&FORMAT=TLE"

# WGS84 Earth ellipsoid
EARTH_A = 6378.137                  # equatorial radius, km
EARTH_F = 1 / 298.257223563         # flattening
EARTH_E2 = EARTH_F * (2 - EARTH_F)  # eccentricity squared

STEP_SECONDS = 5


# ---------------------------------------------------------------- TLE input

def fetch_tle(catnr):
    """Download the latest TLE for a NORAD catalog number from CelesTrak."""
    url = CELESTRAK_URL.format(catnr=catnr)
    with urllib.request.urlopen(url, timeout=20) as response:
        text = response.read().decode()
    return parse_tle(text)


def parse_tle(text):
    """Return (name, line1, line2) from TLE text with or without a name line."""
    lines = [ln.rstrip() for ln in text.splitlines() if ln.strip()]
    line1 = next((ln for ln in lines if ln.startswith("1 ")), None)
    line2 = next((ln for ln in lines if ln.startswith("2 ")), None)
    if line1 is None or line2 is None:
        raise ValueError("No valid TLE found. Check the catalog number or file.")
    name = lines[0].strip() if not lines[0].startswith("1 ") else "Satellite"
    return name, line1, line2


# ------------------------------------------------------------ time and frames

def julian_date(unix_seconds):
    """Split Unix time (UTC seconds) into Julian day number and day fraction."""
    days = np.asarray(unix_seconds, dtype=float) / 86400.0
    whole = np.floor(days)
    return 2440587.5 + whole, days - whole


def gmst(jd, fr):
    """Greenwich Mean Sidereal Time in radians (IAU 1982 model)."""
    t = (jd - 2451545.0 + fr) / 36525.0
    seconds = (67310.54841
               + (876600.0 * 3600.0 + 8640184.812866) * t
               + 0.093104 * t**2
               - 6.2e-6 * t**3)
    return np.radians((seconds % 86400.0) / 240.0)


def teme_to_ecef(r_teme, theta):
    """Rotate TEME positions (n x 3, km) into the Earth-fixed frame."""
    cos_t, sin_t = np.cos(theta), np.sin(theta)
    x = cos_t * r_teme[:, 0] + sin_t * r_teme[:, 1]
    y = -sin_t * r_teme[:, 0] + cos_t * r_teme[:, 1]
    return np.column_stack([x, y, r_teme[:, 2]])


def observer_ecef(lat_deg, lon_deg, alt_km):
    """Earth-fixed position (km) of a point given by geodetic coordinates."""
    lat, lon = np.radians(lat_deg), np.radians(lon_deg)
    n = EARTH_A / np.sqrt(1 - EARTH_E2 * np.sin(lat) ** 2)
    return np.array([
        (n + alt_km) * np.cos(lat) * np.cos(lon),
        (n + alt_km) * np.cos(lat) * np.sin(lon),
        (n * (1 - EARTH_E2) + alt_km) * np.sin(lat),
    ])


def look_angles(r_ecef, lat_deg, lon_deg, alt_km):
    """Azimuth (deg), elevation (deg) and range (km) from observer to satellite."""
    lat, lon = np.radians(lat_deg), np.radians(lon_deg)
    d = r_ecef - observer_ecef(lat_deg, lon_deg, alt_km)

    east = -np.sin(lon) * d[:, 0] + np.cos(lon) * d[:, 1]
    north = (-np.sin(lat) * np.cos(lon) * d[:, 0]
             - np.sin(lat) * np.sin(lon) * d[:, 1]
             + np.cos(lat) * d[:, 2])
    up = (np.cos(lat) * np.cos(lon) * d[:, 0]
          + np.cos(lat) * np.sin(lon) * d[:, 1]
          + np.sin(lat) * d[:, 2])

    rng = np.linalg.norm(d, axis=1)
    az = np.degrees(np.arctan2(east, north)) % 360.0
    el = np.degrees(np.arcsin(up / rng))
    return az, el, rng


# ------------------------------------------------------------- pass finding

def propagate(line1, line2, unix_times):
    """Run SGP4 and return Earth-fixed positions (n x 3, km)."""
    from sgp4.api import Satrec

    sat = Satrec.twoline2rv(line1, line2)
    jd, fr = julian_date(unix_times)
    error, r_teme, _ = sat.sgp4_array(jd, fr)
    if np.any(error != 0):
        raise RuntimeError("SGP4 failed. The TLE may be too old or invalid.")

    epoch_unix = (sat.jdsatepoch - 2440587.5 + sat.jdsatepochF) * 86400.0
    return teme_to_ecef(r_teme, gmst(jd, fr)), epoch_unix


def find_passes(unix_times, az, el, min_el):
    """Find every complete pass where elevation rises above min_el."""
    above = el >= min_el
    edges = np.diff(above.astype(int))
    rises = np.where(edges == 1)[0] + 1
    sets = np.where(edges == -1)[0]

    passes = []
    for start in rises:
        later = sets[sets >= start]
        if len(later) == 0:
            break                      # pass is cut off by the end of the window
        end = later[0]
        peak = start + int(np.argmax(el[start:end + 1]))
        passes.append({
            "rise": unix_times[start], "rise_az": az[start],
            "peak": unix_times[peak], "peak_az": az[peak], "peak_el": el[peak],
            "set": unix_times[end], "set_az": az[end],
            "slice": slice(start, end + 1),
        })
    return passes


def compass(az_deg):
    """Convert an azimuth in degrees to a compass direction."""
    names = ["N", "NE", "E", "SE", "S", "SW", "W", "NW"]
    return names[int((az_deg + 22.5) // 45) % 8]


# ------------------------------------------------------------------- output

def print_passes(passes, name, tz):
    def fmt(unix, pattern):
        return datetime.fromtimestamp(unix, tz).strftime(pattern)

    print(f"\nPasses of {name}\n")
    print(f"{'Date':<12}{'Rise':<10}{'Peak':<10}{'Set':<10}"
          f"{'Max el':>7}{'Duration':>10}   Direction")
    print("-" * 72)
    for p in passes:
        minutes, seconds = divmod(int(p["set"] - p["rise"]), 60)
        path = f"{compass(p['rise_az'])} to {compass(p['set_az'])}"
        print(f"{fmt(p['rise'], '%a %d %b'):<12}"
              f"{fmt(p['rise'], '%H:%M:%S'):<10}"
              f"{fmt(p['peak'], '%H:%M:%S'):<10}"
              f"{fmt(p['set'], '%H:%M:%S'):<10}"
              f"{p['peak_el']:>6.1f}°"
              f"{minutes:>6d}m {seconds:02d}s   {path}")
    if not passes:
        print("No passes found in this window.")


def plot_pass(p, az, el, name, tz, path):
    """Save a sky chart of one pass: center is straight up, edge is the horizon."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=(6, 6.4), dpi=150)
    ax = fig.add_subplot(projection="polar")
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)

    ax.plot(np.radians(az[p["slice"]]), 90 - el[p["slice"]], linewidth=2.5)
    ax.plot(np.radians(p["rise_az"]), 90 - el[p["slice"]][0], "o", label="Rise")
    ax.plot(np.radians(p["peak_az"]), 90 - p["peak_el"], "*", markersize=14,
            label=f"Peak {p['peak_el']:.0f}°")
    ax.plot(np.radians(p["set_az"]), 90 - el[p["slice"]][-1], "s", label="Set")

    ax.set_rlim(0, 90)
    ax.set_yticks([0, 30, 60, 90])
    ax.set_yticklabels(["90°", "60°", "30°", "0°"])
    ax.set_xticks(np.radians(np.arange(0, 360, 45)))
    ax.set_xticklabels(["N", "NE", "E", "SE", "S", "SW", "W", "NW"])

    when = datetime.fromtimestamp(p["rise"], tz).strftime("%a %d %b %Y, %H:%M")
    ax.set_title(f"{name}\n{when}", pad=18)
    ax.legend(loc="lower left", bbox_to_anchor=(-0.12, -0.12))
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


# --------------------------------------------------------------------- main

def main():
    parser = argparse.ArgumentParser(description="Predict satellite passes over a location.")
    parser.add_argument("--catnr", type=int, default=25544,
                        help="NORAD catalog number (default 25544, the ISS)")
    parser.add_argument("--tle-file", help="read the TLE from a file instead of downloading")
    parser.add_argument("--lat", type=float, default=21.2703, help="latitude in degrees (default Taif)")
    parser.add_argument("--lon", type=float, default=40.4158, help="longitude in degrees (default Taif)")
    parser.add_argument("--alt", type=float, default=1879.0, help="altitude in meters (default Taif)")
    parser.add_argument("--days", type=float, default=3.0, help="how many days ahead to search")
    parser.add_argument("--min-el", type=float, default=10.0, help="minimum elevation in degrees")
    parser.add_argument("--utc-offset", type=float, default=3.0, help="local time offset from UTC in hours")
    parser.add_argument("--start", help='search start in UTC, "YYYY-MM-DD HH:MM" (default now)')
    parser.add_argument("--plot", default="best_pass.png", help="where to save the sky chart")
    parser.add_argument("--no-plot", action="store_true", help="skip the sky chart")
    args = parser.parse_args()

    if args.tle_file:
        name, line1, line2 = parse_tle(Path(args.tle_file).read_text())
    else:
        name, line1, line2 = fetch_tle(args.catnr)

    if args.start:
        start = datetime.strptime(args.start, "%Y-%m-%d %H:%M").replace(tzinfo=timezone.utc)
    else:
        start = datetime.now(timezone.utc)

    unix_times = start.timestamp() + np.arange(0, args.days * 86400, STEP_SECONDS)
    r_ecef, epoch_unix = propagate(line1, line2, unix_times)

    tle_age_days = abs(start.timestamp() - epoch_unix) / 86400
    if tle_age_days > 14:
        print(f"Warning: this TLE is {tle_age_days:.0f} days from the search start. "
              "Predictions lose accuracy as a TLE ages.")

    az, el, _ = look_angles(r_ecef, args.lat, args.lon, args.alt / 1000.0)
    passes = find_passes(unix_times, az, el, args.min_el)

    tz = timezone(timedelta(hours=args.utc_offset))
    print_passes(passes, name, tz)

    if passes and not args.no_plot:
        best = max(passes, key=lambda p: p["peak_el"])
        plot_pass(best, az, el, name, tz, args.plot)
        print(f"\nSky chart of the highest pass saved to {args.plot}")


if __name__ == "__main__":
    main()
