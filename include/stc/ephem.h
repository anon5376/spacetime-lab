#ifndef STC_EPHEM_H
#define STC_EPHEM_H

/* SPACETIME LAB — ephem: real Solar System catalog + calendar time for
 * the space simulator.
 *
 * Data (claim class EXTERNAL): barycentric state vectors in the J2000
 * ecliptic frame (ICRF axes rotated to the ecliptic of J2000, as JPL
 * Horizons REF_PLANE=ECLIPTIC) at each catalog epoch, GM values, and
 * mean radii, all taken from NASA JPL Horizons
 * (https://ssd.jpl.nasa.gov/horizons/), fetched once by
 * tools/fetch_horizons.sh into data/horizons/ (raw responses kept,
 * query URLs and fetch date recorded) and compiled into
 * src/space/ephem_data.c by that tool. Nothing is typed in by hand.
 *
 * Planet vs system barycentre: when moons are loaded, each planet is
 * the planet body (e.g. Horizons 599 Jupiter, GM of the planet alone)
 * plus its moons; when moons are NOT loaded, the planet is replaced by
 * its system barycentre (Horizons 5, GM of planet + moons) so the
 * outer dynamics stay right. Earth: 399 + 301 Moon, or 3 (Earth-Moon
 * barycentre).
 *
 * Time: TDB throughout. t = seconds since J2000.0 = JD 2451545.0 TDB.
 * Calendar conversion is proleptic Gregorian (valid for any JD > 0);
 * the TDB-UTC difference (~69 s today) is NOT applied: dates shown are
 * TDB calendar dates.
 */

#include "stc/space.h"

/* Catalog groups (bitmask). */
#define STC_EPH_CORE 1u      /* Sun + 8 planets (barycentres if no moons) */
#define STC_EPH_MOONS 2u     /* major moons (Moon, Galileans, Titan, ...) */
#define STC_EPH_DWARFS 4u    /* Pluto(+Charon w/ moons), Eris, Makemake,
                                Haumea, Ceres */
#define STC_EPH_ASTEROIDS 8u /* largest main-belt asteroids (Vesta,
                                Pallas, Hygiea, ...) as massive bodies */
#define STC_EPH_ALL 15u

typedef struct {
    const char *name;     /* display name, e.g. "Jupiter" */
    int horizons_id;      /* Horizons COMMAND id used (e.g. 599) */
    int system_id;        /* barycentre id used when moons are off
                             (e.g. 5), or 0 if not applicable */
    int kind;             /* STC_SK_* */
    const char *parent;   /* primary's catalog name, NULL for the Sun */
    unsigned group;       /* STC_EPH_* */
    double gm;            /* m^3/s^2 of the body itself */
    double gm_system;     /* m^3/s^2 of the system barycentre (or gm) */
    double radius;        /* mean radius, m */
    /* Zonal gravity (EXTERNAL, cited in data/horizons/MANIFEST.md);
     * all 0 when not modelled. Applied by stc_ephem_load when MOONS are
     * loaded (the planet is then a separate body from its moons). */
    double j2;
    double j4;
    double rref;          /* reference radius for J2/J4, m */
    double pole[3];       /* spin pole unit vector, J2000 ecliptic */
} StcEphemBody;

/* Number of catalog entries, entry i (NULL out of range), lookup. */
int stc_ephem_count(void);
const StcEphemBody *stc_ephem_get(int i);
int stc_ephem_find(const char *name);

/* Catalog epochs (JD TDB) with full state vectors for every entry.
 * Epoch 0 is J2000.0; others are listed in ephem_data.c. */
int stc_ephem_epoch_count(void);
double stc_ephem_epoch_jd(int k);

/* Horizons state of catalog entry i at epoch k: barycentric J2000
 * ecliptic, m and m/s. use_system != 0 returns the system-barycentre
 * vector (for entries with system_id). Returns 0 or -1. */
int stc_ephem_state(int i, int k, int use_system, double pos[3],
                    double vel[3]);

/* Provenance line: data source, DE version reported by Horizons, fetch
 * date. Never NULL. */
const char *stc_ephem_source(void);

/* Load the groups at epoch k into s (cleared first; s->t set to the
 * epoch). Applies the planet/barycentre rule above, sets parents,
 * kinds, radii. belt_tracers > 0 adds that many main-belt test
 * particles on osculating orbits about the Sun drawn from a seeded
 * distribution (a 2.1-3.3 AU, e < 0.3, i < 20 deg, uniform angles):
 * ILLUSTRATIVE (a sampled population, not real asteroids).
 * Returns the number of massive bodies loaded, -1 on bad args,
 * -2 on capacity (reserve first). */
int stc_ephem_load(StcSpace *s, int k, unsigned groups, int belt_tracers,
                   unsigned long long seed);

/* --- time --------------------------------------------------------- */

#define STC_JD_J2000 2451545.0

double stc_ephem_jd(double t);          /* JD TDB from t */
double stc_ephem_t(double jd);          /* t from JD TDB */

/* JD -> proleptic Gregorian calendar (TDB). */
void stc_ephem_cal(double jd, int *year, int *month, int *day, int *hour,
                   int *minute, double *second);

/* Calendar -> JD. Accepts out-of-range fields (normalised). */
double stc_ephem_jd_from_cal(int year, int month, int day, int hour,
                             int minute, double second);

/* Format "YYYY-MM-DD HH:MM:SS TDB" (negative years as -YYYY) into buf.
 * Returns buf. */
char *stc_ephem_format(double t, char *buf, int n);

#endif /* STC_EPHEM_H */
