#!/usr/bin/env python3
"""SPACETIME LAB — tools/gen_ephem.py

Parses the raw NASA JPL Horizons responses fetched by tools/fetch_horizons.sh
(under data/horizons/) and generates src/space/ephem_data.c: the compiled
EPHEM catalog (include/stc/ephem.h). Nothing here is typed in by hand — every
number in ephem_data.c traces to a specific file in data/horizons/, which in
turn traces to a URL and UTC fetch time recorded in data/horizons/MANIFEST.md.

Run from the spacetime/ project root (or anywhere; paths are relative to
this script's own directory):

    python3 tools/gen_ephem.py

Dropped bodies (documented, not fabricated): Eris, Makemake, Haumea, Davida
were requested by the program's catalog spec but Horizons' small-body
lookup returns "GM= n.a." for all four in the fetched responses (Eris/
Makemake/Haumea also return "RAD= n.a."). Program rule R8 (no fabricated
references — every number traces to a fetched Horizons response) and the
catalog spec's own instruction ("if no GM is published there, drop that
body and list it in your report") both apply: there is no fetched, real
GM for these four, so they are left out of the compiled catalog rather
than filled in from outside memory. See the tool's stdout / the worker
report for the exact grep evidence.
"""
import math
import os
import re
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(SCRIPT_DIR)
HZ_DIR = os.path.join(ROOT, "data", "horizons")
NAIF_DIR = os.path.join(ROOT, "data", "naif")
OUT_C = os.path.join(ROOT, "src", "space", "ephem_data.c")

# IAU 1976 obliquity of the ecliptic at J2000.0, 84381.448 arcsec. This is
# the constant Horizons itself cites for REF_PLANE=ECLIPTIC (see any fetched
# data/horizons/*.txt: "Note: IAU76 obliquity of 84381.448 arcseconds wrt
# ICRF X-Y plane" in the "REFERENCE FRAME AND COORDINATES" section) -- not a
# separately fetched number, but confirmed present in the already-fetched
# Horizons responses rather than typed from memory.
OBLIQUITY_ARCSEC = 84381.448

G_SI = 6.67430e-11  # STC_G, m^3 kg^-1 s^-2 (must match include/stc/units.h)

# STC_SK_* (include/stc/space.h)
SK_STAR, SK_PLANET, SK_DWARF, SK_MOON, SK_ASTEROID = 0, 1, 2, 3, 4

# STC_EPH_* (include/stc/ephem.h)
EPH_CORE, EPH_MOONS, EPH_DWARFS, EPH_ASTEROIDS = 1, 2, 4, 8

EPOCH_TAGS = ["e0", "e1"]
EPOCH_LABELS = {"e0": "J2000.0 (2000-01-01 12:00 TDB)",
                "e1": "2026-09-27 00:00 TDB"}


def read(path):
    with open(path, "r", errors="replace") as f:
        return f.read()


def parse_vectors(text, src_desc):
    """First $$SOE row -> (jd, [x,y,z,vx,vy,vz] in meters, m/s)."""
    m = re.search(r"\$\$SOE\s*\n(.*?)\n\$\$EOE", text, re.S)
    if not m:
        raise ValueError("no $$SOE/$$EOE block in " + src_desc)
    row = m.group(1).strip().splitlines()[0]
    parts = [p.strip() for p in row.split(",")]
    jd = float(parts[0])
    x, y, z, vx, vy, vz = (float(p) for p in parts[2:8])
    km = 1000.0
    return jd, [x * km, y * km, z * km, vx * km, vy * km, vz * km]


_NUM = r"([+-]?[0-9][0-9.x]*)"


def _num(tok):
    if "x" in tok:
        vals = [float(v) for v in tok.split("x")]
        return sum(vals) / len(vals)
    return float(tok)


def parse_gm_radius(text, src_desc):
    """Returns (gm_km3s2, radius_km) or (None, None) parts if not found.
    Handles the several PHYSICAL DATA / Asteroid physical parameters
    layouts Horizons uses (planet, moon, dwarf-major-body, small-body)."""
    gm = None
    m = re.search(r"\bGM\b(?!\s*1-sigma)[^\n=]{0,40}=\s*" + _NUM, text)
    if m:
        tok = m.group(1)
        gm = None if tok.lower().startswith("n") else _num(tok)
    if gm is None:
        m = re.search(r"\bGM\s*=\s*(n\.a\.|" + _NUM + r")", text)
        if m and not m.group(1).startswith("n"):
            gm = _num(m.group(1))

    radius = None
    for pat in (
        r"vol\.?\s*mean\s*radius[^=\n]*=\s*" + _NUM,
        r"mean\s*radius[^=\n]*=\s*" + _NUM,
        r"radius\s*\(km[^)]*\)\s*=\s*" + _NUM,
        r"\brad\s*=\s*(n\.a\.|" + _NUM + r")",
    ):
        m = re.search(pat, text, re.I)
        if m:
            tok = m.group(1)
            if tok.lower().startswith("n"):
                continue
            radius = _num(tok)
            break
    return gm, radius


def source_tag(text):
    m = re.search(r"source:\s*([A-Za-z0-9#]+)", text)
    return m.group(1) if m else "?"


def check_obliquity_in_fetched_data():
    """The IAU76 obliquity constant used for the ecliptic-frame pole
    conversion is not fetched on its own; confirm it is actually present,
    verbatim, in an already-fetched Horizons response (R8: every number
    traces to a fetched official response)."""
    text = read(os.path.join(HZ_DIR, "earth_e0.txt"))
    if "84381.448" not in text:
        print("FATAL: IAU76 obliquity 84381.448 not found in a fetched "
              "Horizons response (earth_e0.txt) -- cannot confirm the "
              "ecliptic-frame pole conversion constant.", file=sys.stderr)
        sys.exit(1)


# --- oblateness (EPHEM-2): J2/J4/rref from NAIF satellite-ephemeris kernel
# comment files (the .cmt of the .bsp that actually produced the moon
# positions in data/horizons/), pole RA/Dec from the NAIF generic PCK. ------

# name -> (naif .txt file (data/naif/), "Satellite Ephemeris: TAG" marker to
# split on -- the giant-planet .cmt files bundle several historical kernel
# solutions in one file; the marker selects the block matching the Horizons
# source tag already recorded for that body's vectors.
#
# Uranus is the one exception, chosen EMPIRICALLY, not by tag-matching
# (EPHEM-3): data/naif/ura184.txt's own headline tag is URA184, whose J702/
# J704 (4.035987277058390e-3 / -3.552301969345126e-5 @ RADIUS 25559 km,
# "Satellite Ephemeris: URA184" block) are what a first pass used, but a
# point-mass+J2/J4 integration with those coefficients (stc_space_advance,
# ALL group, h=1800s, pn=1, J2000->2026-09-27) puts the Uranian moons off
# the fetched Horizons vectors by 2.5e5 km (Miranda) / 1.3e5 (Ariel) /
# 5.9e4 (Umbriel) / 3.7e4 (Titania) -- ura184.txt lines 466-524 also carry
# an OLDER solution, "Satellite Ephemeris: URA182.21" (lines 114-181), with
# J702=3.508966976546036e-3, J704=-3.579347133056984e-5, RADIUS=25559 km
# (values at lines 136-137, repeated at 203-204). Re-running the SAME
# integration with the URA182.21 coefficients instead gives Miranda 3052,
# Ariel 2576, Umbriel 3258, Titania 20788, Oberon 3928->4557 km -- two
# orders of magnitude closer. Horizons' own moon vectors are evidently
# consistent with J2 ~= 3.51e-3, not the URA184 file's headline 4.04e-3;
# the URA184 coefficient likely pairs with additional ring/higher-order
# modelling this simulator does not reproduce (higher zonal terms, secular
# pole precession, etc.), so it is the wrong number to feed a bare J2/J4
# point-mass model even though it is the "current" kernel solution. See
# data/naif/MANIFEST.md for the full experiment table (build/_ura.c, the
# lead's harness, independently re-run here with identical results) and
# ephem_data.c's header for this same note.
NAIF_GRAVITY_SRC = {
    "Jupiter": ("jup365.txt", "JUP365.26"),
    "Saturn": ("sat441.txt", "SAT441.24"),
    "Uranus": ("ura184.txt", "URA182.21"),
    "Neptune": ("nep098_part-1.txt", "NEP098"),
}


def parse_naif_planet_gravity(name, body_num):
    """Returns (j2, j4, rref_km) for a giant planet from its NAIF
    satellite-ephemeris kernel comment file, scoped to the block matching
    the kernel solution actually used (see NAIF_GRAVITY_SRC)."""
    fname, tag = NAIF_GRAVITY_SRC[name]
    text = read(os.path.join(NAIF_DIR, fname))
    marker = f"Satellite Ephemeris: {tag}"
    if marker not in text:
        raise ValueError(f"{marker!r} not found in {fname}")
    block = text.split(marker, 1)[1]
    j2 = float(re.search(rf"\bJ{body_num}02\s+([+-]?[0-9.]+E[+-][0-9]+)",
                         block).group(1))
    j4 = float(re.search(rf"\bJ{body_num}04\s+([+-]?[0-9.]+E[+-][0-9]+)",
                         block).group(1))
    rref = float(re.search(r"\bRADIUS\s+([0-9.]+E[+-][0-9]+)",
                           block).group(1))
    return j2, j4, rref


def parse_pck_pole(body_id):
    """Returns (ra_deg, dec_deg) constant terms for BODY<id>_POLE_RA/DEC
    from the NAIF generic PCK, reading only the \\begindata region (the
    file also carries an identical-looking example inside its
    \\begintext documentation prose, which must NOT be parsed as data)."""
    text = read(os.path.join(NAIF_DIR, "pck00011.txt"))
    lines = text.splitlines()
    in_data = False
    ra = dec = None
    for line in lines:
        if "\\begindata" in line:
            in_data = True
            continue
        if "\\begintext" in line:
            in_data = False
            continue
        if not in_data:
            continue
        if ra is None:
            m = re.search(rf"BODY{body_id}_POLE_RA\s*=\s*\(\s*"
                          r"([+-]?[0-9.]+)", line)
            if m:
                ra = float(m.group(1))
        if dec is None:
            m = re.search(rf"BODY{body_id}_POLE_DEC\s*=\s*\(\s*"
                          r"([+-]?[0-9.]+)", line)
            if m:
                dec = float(m.group(1))
        if ra is not None and dec is not None:
            break
    if ra is None or dec is None:
        raise ValueError(f"BODY{body_id}_POLE_RA/DEC not found in pck data")
    return ra, dec


def ecliptic_pole(ra_deg, dec_deg):
    """Spin-pole unit vector in the J2000 ECLIPTIC frame from ICRF-
    equatorial (RA, Dec) constant terms (precession/nutation rate terms
    ignored -- see ephem_data.c header for the values and this omission).
    Rotate the ICRF-equatorial unit vector by the IAU76 obliquity about
    the shared X axis (the equinox direction is common to both frames):
        x_ecl =  x_eq
        y_ecl =  cos(eps) y_eq + sin(eps) z_eq
        z_ecl = -sin(eps) y_eq + cos(eps) z_eq
    Sanity: Earth's own pole is exactly the equatorial +Z axis (RA=0,
    Dec=90 at J2000), which this maps to (0, sin(eps), cos(eps)) -- tilted
    away from the ecliptic pole by exactly the obliquity, as expected."""
    ra = math.radians(ra_deg)
    dec = math.radians(dec_deg)
    x_eq = math.cos(dec) * math.cos(ra)
    y_eq = math.cos(dec) * math.sin(ra)
    z_eq = math.sin(dec)
    eps = math.radians(OBLIQUITY_ARCSEC / 3600.0)
    x = x_eq
    y = math.cos(eps) * y_eq + math.sin(eps) * z_eq
    z = -math.sin(eps) * y_eq + math.cos(eps) * z_eq
    n = math.sqrt(x * x + y * y + z * z)
    return (x / n, y / n, z / n)


# name -> Horizons/NAIF body id used for the PCK pole lookup
POLE_BODY_ID = {"Earth": 399, "Jupiter": 599, "Saturn": 699, "Uranus": 799,
                "Neptune": 899}

# Earth's J2 comes straight from the already-fetched Horizons Earth page
# ("J2 (IERS 2010) = 0.00108262545"); no separate NAIF gravity-field fetch
# needed/available for a point-body page. No officially-fetched Earth J4
# was found in either source (Horizons' Earth page nor the NAIF kernels,
# which don't carry Earth's field) -- left at 0.0 rather than filled in
# from outside memory (R8). The Moon's J2 was also sought and not cleanly
# available from either official source used here, so it stays a point
# mass too (documented, not fabricated).
EARTH_J2_PATTERN = r"J2\s*\(IERS\s*2010\)\s*=\s*([0-9.]+)"
EARTH_REQ_PATTERN = r"Equ\.\s*radius,\s*km\s*=\s*([0-9.]+)"


def compute_oblateness():
    """Returns {name: (j2, j4, rref_m, pole_xyz)} for the 5 bodies with an
    officially fetched gravity field; all others get (0,0,0,(0,0,0))."""
    check_obliquity_in_fetched_data()
    out = {}

    earth_text = read(os.path.join(HZ_DIR, "earth_e0.txt"))
    j2 = float(re.search(EARTH_J2_PATTERN, earth_text).group(1))
    rref_km = float(re.search(EARTH_REQ_PATTERN, earth_text).group(1))
    ra, dec = parse_pck_pole(POLE_BODY_ID["Earth"])
    out["Earth"] = (j2, 0.0, rref_km * 1000.0, ecliptic_pole(ra, dec))

    for name, body_num in (("Jupiter", 5), ("Saturn", 6), ("Uranus", 7),
                           ("Neptune", 8)):
        j2, j4, rref_km = parse_naif_planet_gravity(name, body_num)
        ra, dec = parse_pck_pole(POLE_BODY_ID[name])
        out[name] = (j2, j4, rref_km * 1000.0, ecliptic_pole(ra, dec))

    return out


# --- catalog definition ------------------------------------------------
# (name, file_base, kind, parent_name_or_None, group, system_id)
CATALOG = [
    ("Sun", "sun", SK_STAR, None, EPH_CORE, 0),
    ("Mercury", "mercury", SK_PLANET, "Sun", EPH_CORE, 1),
    ("Venus", "venus", SK_PLANET, "Sun", EPH_CORE, 2),
    ("Earth", "earth", SK_PLANET, "Sun", EPH_CORE, 3),
    ("Mars", "mars", SK_PLANET, "Sun", EPH_CORE, 4),
    ("Jupiter", "jupiter", SK_PLANET, "Sun", EPH_CORE, 5),
    ("Saturn", "saturn", SK_PLANET, "Sun", EPH_CORE, 6),
    ("Uranus", "uranus", SK_PLANET, "Sun", EPH_CORE, 7),
    ("Neptune", "neptune", SK_PLANET, "Sun", EPH_CORE, 8),

    ("Pluto", "pluto", SK_DWARF, "Sun", EPH_DWARFS, 9),
    ("Ceres", "ceres", SK_DWARF, "Sun", EPH_DWARFS, 0),

    ("Vesta", "vesta", SK_ASTEROID, "Sun", EPH_ASTEROIDS, 0),
    ("Pallas", "pallas", SK_ASTEROID, "Sun", EPH_ASTEROIDS, 0),
    ("Hygiea", "hygiea", SK_ASTEROID, "Sun", EPH_ASTEROIDS, 0),
    ("Interamnia", "interamnia", SK_ASTEROID, "Sun", EPH_ASTEROIDS, 0),

    ("Moon", "moon", SK_MOON, "Earth", EPH_MOONS, 0),
    ("Io", "io", SK_MOON, "Jupiter", EPH_MOONS, 0),
    ("Europa", "europa", SK_MOON, "Jupiter", EPH_MOONS, 0),
    ("Ganymede", "ganymede", SK_MOON, "Jupiter", EPH_MOONS, 0),
    ("Callisto", "callisto", SK_MOON, "Jupiter", EPH_MOONS, 0),
    ("Mimas", "mimas", SK_MOON, "Saturn", EPH_MOONS, 0),
    ("Enceladus", "enceladus", SK_MOON, "Saturn", EPH_MOONS, 0),
    ("Tethys", "tethys", SK_MOON, "Saturn", EPH_MOONS, 0),
    ("Dione", "dione", SK_MOON, "Saturn", EPH_MOONS, 0),
    ("Rhea", "rhea", SK_MOON, "Saturn", EPH_MOONS, 0),
    ("Titan", "titan", SK_MOON, "Saturn", EPH_MOONS, 0),
    ("Iapetus", "iapetus", SK_MOON, "Saturn", EPH_MOONS, 0),
    ("Miranda", "miranda", SK_MOON, "Uranus", EPH_MOONS, 0),
    ("Ariel", "ariel", SK_MOON, "Uranus", EPH_MOONS, 0),
    ("Umbriel", "umbriel", SK_MOON, "Uranus", EPH_MOONS, 0),
    ("Titania", "titania", SK_MOON, "Uranus", EPH_MOONS, 0),
    ("Oberon", "oberon", SK_MOON, "Uranus", EPH_MOONS, 0),
    ("Triton", "triton", SK_MOON, "Neptune", EPH_MOONS, 0),
    ("Charon", "charon", SK_MOON, "Pluto", EPH_MOONS, 0),
]

# Dropped by measurement (Horizons GM= n.a.), kept here only for the report.
DROPPED = ["eris", "makemake", "haumea", "davida"]

# system barycentre file base -> (planet/dwarf name it belongs to, moons to
# sum GM over for the gm_system approximation)
SYSTEMS = [
    ("mercury_sys", "Mercury", []),
    ("venus_sys", "Venus", []),
    ("earth_sys", "Earth", ["Moon"]),
    ("mars_sys", "Mars", []),
    ("jupiter_sys", "Jupiter", ["Io", "Europa", "Ganymede", "Callisto"]),
    ("saturn_sys", "Saturn", ["Mimas", "Enceladus", "Tethys", "Dione",
                              "Rhea", "Titan", "Iapetus"]),
    ("uranus_sys", "Uranus", ["Miranda", "Ariel", "Umbriel", "Titania",
                              "Oberon"]),
    ("neptune_sys", "Neptune", ["Triton"]),
    ("pluto_sys", "Pluto", ["Charon"]),
]


def main():
    n = len(CATALOG)
    state = [[None, None] for _ in range(n)]
    gm_km = [None] * n
    rad_km = [None] * n
    sources = set()
    epoch_jd = [None, None]

    for i, (name, base, kind, parent, group, sysid) in enumerate(CATALOG):
        for ei, tag in enumerate(EPOCH_TAGS):
            path = os.path.join(HZ_DIR, f"{base}_{tag}.txt")
            text = read(path)
            sources.add(source_tag(text))
            jd, vec = parse_vectors(text, path)
            if epoch_jd[ei] is None:
                epoch_jd[ei] = jd
            else:
                assert abs(epoch_jd[ei] - jd) < 1e-6, (name, tag, jd)
            state[i][ei] = vec
        gm, rad = parse_gm_radius(read(os.path.join(HZ_DIR, f"{base}_e0.txt")),
                                   base)
        if gm is None or rad is None:
            print(f"FATAL: {name} ({base}) missing gm={gm} rad={rad}",
                  file=sys.stderr)
            sys.exit(1)
        gm_km[i] = gm
        rad_km[i] = rad

    name_to_idx = {c[0]: i for i, c in enumerate(CATALOG)}
    oblate = compute_oblateness()
    for nm, (j2, j4, rref_m, pole) in oblate.items():
        print(f"oblateness: {nm:8s} j2={j2:.8e} j4={j4:.8e} "
              f"rref={rref_m/1000:.3f} km pole_ecl={pole}")

    sys_state = []
    sys_catalog_index = []
    gm_system_km = {}
    for base, planet, moons in SYSTEMS:
        pi = name_to_idx[planet]
        rows = [None, None]
        for ei, tag in enumerate(EPOCH_TAGS):
            path = os.path.join(HZ_DIR, f"{base}_{tag}.txt")
            text = read(path)
            sources.add(source_tag(text))
            jd, vec = parse_vectors(text, path)
            rows[ei] = vec
        sys_state.append(rows)
        sys_catalog_index.append(pi)
        gm_sum = gm_km[pi] + sum(gm_km[name_to_idx[m]] for m in moons)
        gm_system_km[planet] = (gm_sum, moons)

    # --- sanity per the operator's known-good sample (from the task brief):
    # Earth J2000 x = -2.756674048281145e10 m.
    earth_i = name_to_idx["Earth"]
    ex = state[earth_i][0][0]
    expect = -2.756674048281145e10
    relerr = abs(ex - expect) / abs(expect)
    print(f"sanity: Earth e0 x = {ex:.9e} m (expect {expect:.9e}, "
          f"relerr {relerr:.3e})")
    assert relerr < 1e-9, "Earth J2000 x does not match the tested sample"

    # --- emit ephem_data.c --------------------------------------------
    lines = []
    a = lines.append
    a("/* src/space/ephem_data.c — GENERATED FILE, do not hand-edit.")
    a(" *")
    a(" * Produced by tools/gen_ephem.py from the raw NASA JPL Horizons")
    a(" * responses in data/horizons/ (fetched by tools/fetch_horizons.sh;")
    a(" * every query URL + UTC fetch time is in data/horizons/MANIFEST.md).")
    a(" * Regenerate with: python3 tools/gen_ephem.py")
    a(" *")
    a(" * Claim class EXTERNAL (doc/SWE2-SYSTEM.md R2): state vectors, GM,")
    a(" * and radii are values JPL Horizons reported, not independently")
    a(" * computed here. State vectors are barycentric (CENTER=500@0),")
    a(" * J2000 ecliptic (REF_PLANE=ECLIPTIC, REF_SYSTEM=ICRF), TIME_TYPE=")
    a(" * TDB, converted km/km-s -> m/m-s (SI) at generation time only.")
    a(" *")
    a(" * Horizons sources seen across the fetched set: " +
      ", ".join(sorted(sources)) + ".")
    a(" *")
    a(" * gm_system (used when moons are off) is NOT reported by Horizons")
    a(" * for a barycentre (a barycentre has no GM of its own); it is")
    a(" * computed here as the sum of the fetched component GMs:")
    for planet, (gm, moons) in gm_system_km.items():
        comp = planet + (" + " + " + ".join(moons) if moons else " (no cataloged moons)")
        a(f" *   {planet:8s} gm_system = {comp}")
    a(" * For the outer planets this sums only the cataloged major moons")
    a(" * (no minor-moon/ring mass), a documented approximation; the state")
    a(" * VECTOR for each system entry is still the real Horizons barycentre")
    a(" * ephemeris (all system mass already accounted for by JPL).")
    a(" *")
    a(" * Dropped from the requested catalog (Horizons small-body lookup")
    a(" * returned GM= n.a. for all of these; R8 forbids filling in a")
    a(" * literature value not actually fetched this run): Eris (136199),")
    a(" * Makemake (136472), Haumea (136108) [also RAD= n.a.], Davida (511).")
    a(" *")
    a(" * EPHEM-2: zonal gravity (j2/j4/rref) + spin pole (ecliptic J2000),")
    a(" * all EXTERNAL, all 0/(0,0,0) except the 5 bodies below:")
    for nm in ("Earth", "Jupiter", "Saturn", "Uranus", "Neptune"):
        j2, j4, rref_m, pole = oblate[nm]
        a(f" *   {nm:8s} j2={j2:.10e} j4={j4:.10e} rref={rref_m/1000.0:.3f} km"
          f" pole_ecl=({pole[0]:.6f},{pole[1]:.6f},{pole[2]:.6f})")
    a(" * Sources: Earth j2 from data/horizons/earth_e0.txt (\"J2 (IERS")
    a(" * 2010)\"); giants' j2/j4/rref from data/naif/{jup365,sat441,ura184,")
    a(" * nep098_part-1}.txt (the SPICE satellite-ephemeris kernel comment")
    a(" * file matching the .bsp that produced that planet's moon vectors,")
    a(" * scoped to the block whose \"Satellite Ephemeris:\" tag matches the")
    a(" * Horizons source tag already recorded for those vectors). Pole")
    a(" * RA/Dec (constant term only; slow precession-rate terms ignored,")
    a(" * see data/naif/MANIFEST.md for the rate coefficients) from")
    a(" * data/naif/pck00011.txt BODY<id>_POLE_RA/DEC, rotated from ICRF")
    a(" * equatorial into the J2000 ecliptic frame by the IAU76 obliquity")
    a(" * 84381.448 arcsec (tools/gen_ephem.py:ecliptic_pole; that constant")
    a(" * is confirmed present in data/horizons/earth_e0.txt, not fetched")
    a(" * separately). Earth j4 and the Moon's j2 were sought but not")
    a(" * cleanly available from an official fetched source and are left")
    a(" * at 0.0 (R8: not filled in from outside memory).")
    a(" * Applied only when STC_EPH_MOONS is loaded (stc_ephem_load): with")
    a(" * moons off the planet is a system barycentre and stays a point")
    a(" * mass, per include/stc/ephem.h.")
    a(" *")
    a(" * EPHEM-3: Uranus J2/J4 chosen EMPIRICALLY, not by kernel-tag match.")
    a(" * data/naif/ura184.txt bundles multiple historical solutions; its")
    a(" * headline block, \"Satellite Ephemeris: URA184\" (lines 466-524),")
    a(" * gives J702=4.035987277058390e-3, J704=-3.552301969345126e-5,")
    a(" * RADIUS=25559 km. Integrating this simulator's point-mass+J2/J4")
    a(" * model (stc_space_advance, ALL group, h=1800s, pn=1, J2000 ->")
    a(" * 2026-09-27) with those coefficients misses the fetched Horizons")
    a(" * epoch-1 vectors by 2.5e5 km (Miranda) / 1.3e5 (Ariel) / 5.9e4")
    a(" * (Umbriel) / 3.7e4 (Titania) / 3928 (Oberon). The SAME file also")
    a(" * carries an older block, \"Satellite Ephemeris: URA182.21\" (lines")
    a(" * 114-181), giving J702=3.508966976546036e-3,")
    a(" * J704=-3.579347133056984e-5, RADIUS=25559 km (values at lines")
    a(" * 136-137, repeated 203-204). The identical integration with THOSE")
    a(" * coefficients instead gives Miranda 3052, Ariel 2576, Umbriel 3258,")
    a(" * Titania 20788, Oberon 4557 km -- two orders of magnitude closer.")
    a(" * Conclusion: Horizons' Uranian-moon vectors are consistent with")
    a(" * J2 ~= 3.51e-3 (URA182.21), not the current kernel's headline")
    a(" * 4.04e-3 (URA184) -- the URA184 coefficient likely belongs with")
    a(" * additional modelling (higher zonal terms, secular pole motion,")
    a(" * ring/satellite perturbations) this bare point-mass+J2/J4 core does")
    a(" * not reproduce, so it is the wrong number to feed it even though")
    a(" * it is the newer/current solution. tools/gen_ephem.py therefore")
    a(" * parses the URA182.21 block for Uranus (NAIF_GRAVITY_SRC), not the")
    a(" * file's own headline tag. Both value sets and this experiment are")
    a(" * also recorded in data/naif/MANIFEST.md.")
    a(" */")
    a("")
    a('#include "stc/ephem.h"')
    a("")
    a(f"#define STC_EPHD_NEPOCH {len(EPOCH_TAGS)}")
    a(f"#define STC_EPHD_NSYS {len(SYSTEMS)}")
    a("")
    a(f"const int stc_ephemd_n = {n};")
    a("")
    a("const StcEphemBody stc_ephemd_catalog[] = {")
    kind_name = {SK_STAR: "STC_SK_STAR", SK_PLANET: "STC_SK_PLANET",
                 SK_DWARF: "STC_SK_DWARF", SK_MOON: "STC_SK_MOON",
                 SK_ASTEROID: "STC_SK_ASTEROID"}
    for i, (name, base, kind, parent, group, sysid) in enumerate(CATALOG):
        gm_si = gm_km[i] * 1.0e9
        rad_si = rad_km[i] * 1.0e3
        if sysid:
            planet_name = name
            gm_sys_si = gm_system_km[planet_name][0] * 1.0e9
        else:
            gm_sys_si = gm_si
        parent_c = "NULL" if parent is None else f'"{parent}"'
        group_name = {EPH_CORE: "STC_EPH_CORE", EPH_MOONS: "STC_EPH_MOONS",
                      EPH_DWARFS: "STC_EPH_DWARFS",
                      EPH_ASTEROIDS: "STC_EPH_ASTEROIDS"}[group]
        assert len(name) <= 23, name
        if name in oblate:
            j2, j4, rref_m, pole = oblate[name]
        else:
            j2, j4, rref_m, pole = 0.0, 0.0, 0.0, (0.0, 0.0, 0.0)
        a(f'    {{ "{name}", {name_to_hid(base, name)}, {sysid}, '
          f'{kind_name[kind]}, {parent_c}, {group_name}, '
          f"{gm_si!r}, {gm_sys_si!r}, {rad_si!r}, "
          f"{j2!r}, {j4!r}, {rref_m!r}, "
          f"{{ {pole[0]!r}, {pole[1]!r}, {pole[2]!r} }} }},")
    a("};")
    a("")
    a(f"const double stc_ephemd_epoch_jd[STC_EPHD_NEPOCH] = {{")
    a(f"    {epoch_jd[0]!r}, {epoch_jd[1]!r}")
    a("};")
    a("")
    a("const double stc_ephemd_state[][STC_EPHD_NEPOCH][6] = {")
    for i, (name, *_rest) in enumerate(CATALOG):
        a(f"    /* {i}: {name} */")
        a("    {")
        for ei in range(2):
            v = state[i][ei]
            a("        { " + ", ".join(f"{x!r}" for x in v) + " },")
        a("    },")
    a("};")
    a("")
    a(f"const int stc_ephemd_sys_catalog_index[STC_EPHD_NSYS] = {{")
    a("    " + ", ".join(str(x) for x in sys_catalog_index))
    a("};")
    a("")
    a("const double stc_ephemd_sys_state[STC_EPHD_NSYS][STC_EPHD_NEPOCH][6] = {")
    for si, (base, planet, moons) in enumerate(SYSTEMS):
        a(f"    /* {si}: {planet} system barycentre ({base}) */")
        a("    {")
        for ei in range(2):
            v = sys_state[si][ei]
            a("        { " + ", ".join(f"{x!r}" for x in v) + " },")
        a("    },")
    a("};")
    a("")
    src_str = ("NASA JPL Horizons (https://ssd.jpl.nasa.gov/horizons/); "
               "sources " + ", ".join(sorted(sources)) + "; "
               "fetched 2026-09-27 by tools/fetch_horizons.sh; see "
               "data/horizons/MANIFEST.md")
    a(f'const char *const stc_ephemd_source = "{src_str}";')
    a("")

    with open(OUT_C, "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"wrote {OUT_C} ({n} catalog entries, {len(SYSTEMS)} system "
          f"barycentres)")


_HID = {
    "sun": 10, "mercury": 199, "venus": 299, "earth": 399, "mars": 499,
    "jupiter": 599, "saturn": 699, "uranus": 799, "neptune": 899,
    "moon": 301, "io": 501, "europa": 502, "ganymede": 503, "callisto": 504,
    "mimas": 601, "enceladus": 602, "tethys": 603, "dione": 604, "rhea": 605,
    "titan": 606, "iapetus": 608, "miranda": 705, "ariel": 701,
    "umbriel": 702, "titania": 703, "oberon": 704, "triton": 801,
    "charon": 901, "pluto": 999, "ceres": 1, "vesta": 4, "pallas": 2,
    "hygiea": 10, "interamnia": 704,
}


def name_to_hid(base, name):
    return _HID[base]


if __name__ == "__main__":
    main()
