# Tonight's Sky — Planets + Bright Stars + Constellations
# by Myranda Rose Alvarez

# ----- Importing Necessary Impports ------
import argparse
from datetime import datetime, timezone

from geopy.geocoders import Nominatim
from geopy.exc import GeocoderUnavailable, GeocoderServiceError

from skyfield.api import load, wgs84, Star
from skyfield.data import hipparcos

from astropy.coordinates import SkyCoord, get_constellation
import astropy.units as u


# ----- CLI ------
def parse_args():
    p = argparse.ArgumentParser(
        description="Planets + bright stars for a place/time, with constellations."
    )
    where = p.add_mutually_exclusive_group(required=True)
    where.add_argument("--city", help='City name, e.g. "Amarillo, TX" or "Paris, France"')
    where.add_argument("--lat", type=float, help="Latitude in degrees (e.g., 35.1991)")
    p.add_argument("--lon", type=float, help="Longitude in degrees (required with --lat)")

    when = p.add_mutually_exclusive_group(required=True)
    when.add_argument("--when", help='Local time "YYYY-MM-DD HH:MM"')
    when.add_argument("--now", action="store_true", help="Use current local time")

    p.add_argument("--min-alt", type=float, default=5.0, help="Min altitude (deg) to list")
    p.add_argument("--mag-limit", type=float, default=2.5, help="Max star mag (smaller = brighter)")
    p.add_argument("--top", type=int, default=20, help="Max stars to list")
    p.add_argument("--with-moon", action="store_true", help="Include the Moon")
    p.add_argument("--refraction", action="store_true", help="Apply simple refraction to alt/az")
    return p.parse_args()


# ----- Location and Time ------
def geocode_city(name: str):
    geo = Nominatim(user_agent="tonights_sky")
    try:
        loc = geo.geocode(name, timeout=10)
    except (GeocoderUnavailable, GeocoderServiceError) as e:
        raise SystemExit(
            f"Geocoding failed: {e}\n"
            f"• Quick workaround: pass coordinates with --lat/--lon.\n"
            f"• Permanent macOS fix: run 'Install Certificates.command' for your Python, "
            f"or set SSL_CERT_FILE to certifi’s bundle."
        )
    if not loc:
        raise SystemExit(f"Could not geocode city: {name}")
    return loc.latitude, loc.longitude, loc.address


def resolve_location(args):
    if args.city:
        lat, lon, pretty = geocode_city(args.city)
        return lat, lon, pretty
    if args.lon is None:
        raise SystemExit("--lon is required when using --lat")
    return args.lat, args.lon, f"{args.lat:.4f}, {args.lon:.4f}"


def resolve_time(args):
    if args.now:
        dt_local = datetime.now().astimezone()
    else:
        dt_local = datetime.strptime(args.when, "%Y-%m-%d %H:%M").astimezone()
    dt_utc = dt_local.astimezone(timezone.utc)
    return dt_local, dt_utc


# ----- Day or Night Condition ------
def sun_altitude_degrees(lat, lon, dt_utc, refraction: bool) -> float:
    ts = load.timescale()
    t = ts.from_datetime(dt_utc)
    eph = load("de421.bsp")
    earth = eph["earth"]
    observer = earth + wgs84.latlon(lat, lon)  # <-- anchor to Earth
    sun = eph["sun"]
    alt, az, _ = (observer.at(t).observe(sun).apparent()
                  .altaz(temperature_C=10 if refraction else None,
                         pressure_mbar=1013 if refraction else None))
    return alt.degrees


def twilight_message(sun_alt_deg: float) -> str:
    if sun_alt_deg > 0:
        return f"Daylight (Sun {sun_alt_deg:.1f}°)"
    if sun_alt_deg > -6:
        return f"Civil twilight (Sun {sun_alt_deg:.1f}°)"
    if sun_alt_deg > -12:
        return f"Nautical twilight (Sun {sun_alt_deg:.1f}°)"
    if sun_alt_deg > -18:
        return f"Astronomical twilight (Sun {sun_alt_deg:.1f}°)"
    return f"Night sky (Sun {sun_alt_deg:.1f}°)"


# ----- Constellation + Compass ------
def constellation_of(ra_hours: float, dec_deg: float) -> str:
    c = SkyCoord(ra=ra_hours * u.hourangle, dec=dec_deg * u.deg, frame="icrs")
    return get_constellation(c, short_name=False)


def compass_16(az_deg: float) -> str:
    dirs = ["N","NNE","NE","ENE","E","ESE","SE","SSE",
            "S","SSW","SW","WSW","W","WNW","NW","NNW"]
    idx = int((az_deg + 11.25) // 22.5) % 16
    return dirs[idx]


# ----- Planets + Moon ------
def list_planets(lat, lon, dt_utc, min_alt, with_moon, refraction):
    ts = load.timescale()
    t = ts.from_datetime(dt_utc)
    eph = load("de421.bsp")
    earth = eph["earth"]
    observer = earth + wgs84.latlon(lat, lon) 

    targets = {
        "Mercury": "mercury",
        "Venus": "venus",
        "Mars": "mars",
        "Jupiter": "jupiter barycenter",
        "Saturn": "saturn barycenter",
    }
    if with_moon:
        targets["Moon"] = "moon"

    rows = []
    for label, key in targets.items():
        body = eph[key]
        app = observer.at(t).observe(body).apparent()
        alt, az, _ = app.altaz(temperature_C=10 if refraction else None,
                               pressure_mbar=1013 if refraction else None)
        alt_deg = alt.degrees
        if alt_deg < min_alt:
            continue
        az_deg = az.degrees % 360.0
        ra, dec, _ = app.radec()  # hours, degrees
        const = constellation_of(ra.hours, dec.degrees)
        rows.append((label, alt_deg, az_deg, compass_16(az_deg), const))

    rows.sort(key=lambda r: r[1], reverse=True)
    return rows


# ----- Bright Stars ------
COMMON_NAMES = {
    32349: "Sirius", 30438: "Canopus", 113368: "Arcturus", 91262: "Vega",
    37279: "Procyon", 24436: "Aldebaran", 21421: "Rigel", 27989: "Betelgeuse",
    65474: "Altair", 97649: "Deneb", 24608: "Capella", 100453: "Fomalhaut",
}


def list_bright_stars(lat, lon, dt_utc, mag_limit, top, refraction):
    ts = load.timescale()
    t = ts.from_datetime(dt_utc)
    eph = load("de421.bsp")
    earth = eph["earth"]
    observer = earth + wgs84.latlon(lat, lon)  # <-- anchor to Earth

    # Load Hipparcos
    url = hipparcos.URL
    with load.open(url) as f:
        df = hipparcos.load_dataframe(f)

    bright = df[df["magnitude"] <= mag_limit]
    rows = []
    for hip, row in bright.iterrows():
        star = Star.from_dataframe(row)
        app = observer.at(t).observe(star).apparent()
        alt, az, _ = app.altaz(temperature_C=10 if refraction else None,
                               pressure_mbar=1013 if refraction else None)
        if alt.degrees <= 0:
            continue
        ra_hours = row["ra_degrees"] / 15.0 
        dec_deg = row["dec_degrees"]
        const = constellation_of(ra_hours, dec_deg)
        name = COMMON_NAMES.get(hip, f"HIP {hip}")
        rows.append((name, row["magnitude"], alt.degrees, az.degrees % 360.0, const))

    rows.sort(key=lambda r: (r[1], -r[2]))
    return rows[:top]


# ----- Main ------
def main():
    args = parse_args()
    lat, lon, loc_name = resolve_location(args)
    dt_local, dt_utc = resolve_time(args)

    print(f"\nTonight’s Sky")
    print(f"Location: {loc_name}")
    print(f"Time:     {dt_local:%Y-%m-%d %H:%M %Z}  (UTC {dt_utc:%Y-%m-%d %H:%M})")
    print(
        f"Settings: min-alt={args.min_alt}°, mag-limit={args.mag_limit}, top={args.top}, "
        f"moon={'on' if args.with_moon else 'off'}, refr={'on' if args.refraction else 'off'}"
    )

    # Light condition
    sun_alt = sun_altitude_degrees(lat, lon, dt_utc, args.refraction)
    print("Light:    " + twilight_message(sun_alt) + "\n")

    # Planets
    planet_rows = list_planets(lat, lon, dt_utc, args.min_alt, args.with_moon, args.refraction)
    print("Planets:")
    if not planet_rows:
        print("  (none above min altitude)\n")
    else:
        print(f"  {'Body':<10} {'Alt°':>6} {'Az°':>6}  Dir  Constellation")
        print("  " + "-" * 48)
        for label, alt_deg, az_deg, d, const in planet_rows:
            print(f"  {label:<10} {alt_deg:6.1f} {az_deg:6.1f}  {d:<3}  {const}")
        print()

    # Bright stars
    star_rows = list_bright_stars(lat, lon, dt_utc, args.mag_limit, args.top, args.refraction)
    print(f"Bright stars (mag ≤ {args.mag_limit}):")
    if not star_rows:
        print("  (none above horizon at this magnitude limit)\n")
    else:
        print(f"  {'Star':<20} {'Mag':>5} {'Alt°':>7} {'Az°':>7}  Constellation")
        print("  " + "-" * 60)
        for name, mag, alt_deg, az_deg, const in star_rows:
            print(f"  {name:<20} {mag:5.2f} {alt_deg:7.2f} {az_deg:7.2f}  {const}")
        print()


if __name__ == "__main__":
    main()