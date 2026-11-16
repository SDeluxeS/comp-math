#ifndef BODY_CONSTANTS_H
#define BODY_CONSTANTS_H

#include "mercury_gravity.h"

/*
 * Shared astronomical constants for the RK dynamics and Light-Time parts.
 * Gravitational parameters are expressed in AU^3/day^2, matching the units
 * used by dynamics.c. No GM values are read from SPICE kernels.
 */
#define ASTRO_AU_KM 149597870.7
#define ASTRO_DAY_SEC 86400.0

#define BODY_GM_SUN_AU3_DAY2       2.9591220828411956e-04
#define BODY_GM_MERCURY_AU3_DAY2   MERCURY_GM_AU3_DAY2
#define BODY_GM_VENUS_AU3_DAY2     7.2434523326441187e-10
#define BODY_GM_EARTH_AU3_DAY2     8.8876924467071022e-10
#define BODY_GM_MARS_AU3_DAY2      9.5495488297258119e-11
#define BODY_GM_JUPITER_AU3_DAY2   2.8253458252257917e-07
#define BODY_GM_SATURN_AU3_DAY2    8.4597059933762903e-08
#define BODY_GM_MOON_AU3_DAY2      1.0931894624024351e-11
#define BODY_GM_URANUS_AU3_DAY2    1.2920265649682399e-08
#define BODY_GM_NEPTUNE_AU3_DAY2   1.5243573478851939e-08
#define BODY_GM_PLUTO_AU3_DAY2     2.1750964786022513e-12

/* Convert AU^3/day^2 to km^3/s^2 for formulas using km and seconds. */
#define BODY_GM_AU3_DAY2_TO_KM3_SEC2 \
    ((ASTRO_AU_KM * ASTRO_AU_KM * ASTRO_AU_KM) / \
     (ASTRO_DAY_SEC * ASTRO_DAY_SEC))

#endif /* BODY_CONSTANTS_H */
