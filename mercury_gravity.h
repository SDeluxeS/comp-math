#ifndef MERCURY_GRAVITY_H
#define MERCURY_GRAVITY_H

#ifdef __cplusplus
extern "C" {
#endif

/*
 * GM Mercury from the HgM009.sha header, converted from km^3/s^2
 * to AU^3/day^2. Keep BODY_GM[1] equal to this value.
 */
#define MERCURY_GM_AU3_DAY2 4.912499083282356e-11

/* Load HgM009 coefficients and initialize the Mercury gravity model.
 * Returns 1 on success and a negative SHTOOLS error code on failure.
 */
int mercury_gravity_init(const char *coefficient_file);

/*
 * Compute the non-central (spherical-harmonic) Mercury acceleration.
 * Input:  r_j2000_au, spacecraft position relative to Mercury in J2000 [AU].
 *         et_days, SPICE ET in days from J2000.
 * Output: acceleration in J2000 [AU/day^2].
 * Returns 1 on success, 0 if the model has not been initialized.
 */
int mercury_gravity_acceleration(double et_days,
                                 const double r_j2000_au[3],
                                 double a_j2000_au_day2[3]);

#ifdef __cplusplus
}
#endif

#endif
