#include "mercury_gravity.h"
#include "grav.h"

#include <cmath>

namespace {

constexpr double PI = 3.141592653589793238462643383279502884;
constexpr double DEG_TO_RAD = PI / 180.0;
constexpr double AU_KM = 149597870.7;
constexpr double DAY_SEC = 86400.0;

/*
 * HgM009 uses the gravity-frame orientation consistent with HgM008.
 * Pole coordinates, spin rate and libration amplitude use the HgM008
 * solution (Genova et al., 2019), which HgM009 explicitly follows.
 * Pole precession, W0 and the libration phase convention follow
 * Stark et al. (2015). The periodic coefficient ratios follow the
 * Margot (2009) series; the series is scaled to the HgM008 amplitude.
 */
constexpr double ALPHA0_DEG = 281.0082;
constexpr double DELTA0_DEG = 61.4164;
constexpr double ALPHA_RATE_DEG_PER_CENTURY = -0.032808;
constexpr double DELTA_RATE_DEG_PER_CENTURY = -0.0048464;
constexpr double W0_DEG = 329.5469;
constexpr double W_RATE_DEG_PER_DAY = 6.1385054;
constexpr double MERCURY_MEAN_MOTION_DEG_PER_DAY = 4.09233445;
constexpr double LIBRATION_PHASE_EPOCH_DAYS = 42.71182;
constexpr double LIBRATION_AMPLITUDE_ARCSEC = 40.0;
constexpr double MARGOT_PRIMARY_LIBRATION_ARCSEC = 0.00993822 * 3600.0;

ShGrav grav;
bool grav_initialized = false;

void multiply3(const double A[9], const double B[9], double C[9])
{
    double tmp[9];
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            tmp[3 * i + j] = 0.0;
            for (int k = 0; k < 3; ++k)
                tmp[3 * i + j] += A[3 * i + k] * B[3 * k + j];
        }
    }
    for (int i = 0; i < 9; ++i)
        C[i] = tmp[i];
}

void rotation_z_passive(double angle, double R[9])
{
    const double c = std::cos(angle);
    const double s = std::sin(angle);
    const double values[9] = {
         c,  s, 0.0,
        -s,  c, 0.0,
        0.0, 0.0, 1.0
    };
    for (int i = 0; i < 9; ++i)
        R[i] = values[i];
}

void rotation_x_passive(double angle, double R[9])
{
    const double c = std::cos(angle);
    const double s = std::sin(angle);
    const double values[9] = {
        1.0, 0.0, 0.0,
        0.0,  c,  s,
        0.0, -s,  c
    };
    for (int i = 0; i < 9; ++i)
        R[i] = values[i];
}

/*
 * Transformation from J2000/ICRF to Mercury body-fixed coordinates.
 * With standard right-handed passive coordinate-rotation matrices and
 * alpha=right ascension, delta=declination:
 * R_BF<-J2000 = R3(W) R1(90 deg - delta) R3(90 deg + alpha).
 * This construction makes the third row of R equal to the inertial
 * unit vector toward Mercury's north pole. Matrices below are passive.
 */
void mercury_j2000_to_body_fixed(double et_days, double R[9])
{
    const double T = et_days / 36525.0;

    const double alpha_deg =
        ALPHA0_DEG + ALPHA_RATE_DEG_PER_CENTURY * T;
    const double delta_deg =
        DELTA0_DEG + DELTA_RATE_DEG_PER_CENTURY * T;

    /*
     * Stark et al. (2015): M_k = k*n0*(t+t0), with t in days
     * from J2000. The five amplitudes are the Margot (2009) series.
     */
    const double d = et_days;
    const double base_phase = MERCURY_MEAN_MOTION_DEG_PER_DAY
                            * (d + LIBRATION_PHASE_EPOCH_DAYS);
    const double M1 = 1.0 * base_phase * DEG_TO_RAD;
    const double M2 = 2.0 * base_phase * DEG_TO_RAD;
    const double M3 = 3.0 * base_phase * DEG_TO_RAD;
    const double M4 = 4.0 * base_phase * DEG_TO_RAD;
    const double M5 = 5.0 * base_phase * DEG_TO_RAD;

    const double libration_scale =
        LIBRATION_AMPLITUDE_ARCSEC / MARGOT_PRIMARY_LIBRATION_ARCSEC;
    const double W_libration_deg = libration_scale * (
          0.00993822 * std::sin(M1)
        - 0.00104581 * std::sin(M2)
        - 0.00010280 * std::sin(M3)
        - 0.00002364 * std::sin(M4)
        - 0.00000532 * std::sin(M5));

    const double W_deg = W0_DEG + W_RATE_DEG_PER_DAY * d
                       + W_libration_deg;

    double R3_W[9], R1_theta[9], R3_phi[9], tmp[9];
    rotation_z_passive(W_deg * DEG_TO_RAD, R3_W);
    rotation_x_passive((90.0 - delta_deg) * DEG_TO_RAD, R1_theta);
    rotation_z_passive((90.0 + alpha_deg) * DEG_TO_RAD, R3_phi);

    multiply3(R1_theta, R3_phi, tmp);
    multiply3(R3_W, tmp, R);
}

} // namespace

extern "C" int mercury_gravity_init(const char *coefficient_file)
{
    if (coefficient_file == nullptr)
        return -4;

    const double gm_mercury = MERCURY_GM_AU3_DAY2;
    const double radius_mercury = 2440.0 / AU_KM;

    grav = shGravCreate(160, gm_mercury, radius_mercury);
    const int status = shGravLoadCoefficients(grav, coefficient_file);
    grav_initialized = (status > 0);
    return status;
}

extern "C" int mercury_gravity_acceleration(
    double et_days,
    const double r_j2000_au[3],
    double a_j2000_au_day2[3])
{
    if (!grav_initialized || r_j2000_au == nullptr ||
        a_j2000_au_day2 == nullptr)
        return 0;

    double R_body_from_j2000[9];
    mercury_j2000_to_body_fixed(et_days, R_body_from_j2000);

    /* Position vector: J2000 -> Mercury body-fixed. */
    double r_body_fixed[3] = {0.0, 0.0, 0.0};
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j)
            r_body_fixed[i] += R_body_from_j2000[3 * i + j]
                             * r_j2000_au[j];
    }

    /*
     * HgM009.sha starts at degree 1 (degree 0 is absent), and grav.cilm
     * is zero-initialized. Therefore this is the harmonic correction only;
     * the central point-mass term is already present in the EIH model.
     */
    double a_body_fixed[3];
    shGravCalculate(grav,
                    r_body_fixed[0],
                    r_body_fixed[1],
                    r_body_fixed[2],
                    a_body_fixed);

    /* Acceleration vector: Mercury body-fixed -> J2000 (R transpose). */
    for (int j = 0; j < 3; ++j) {
        a_j2000_au_day2[j] = 0.0;
        for (int i = 0; i < 3; ++i)
            a_j2000_au_day2[j] += R_body_from_j2000[3 * i + j]
                                * a_body_fixed[i];
    }

    return 1;
}
