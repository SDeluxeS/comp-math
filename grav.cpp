#include <cmath>
#include <string.h>
#include "grav.h"
#include "include/shtools.h"

ShGrav shGravCreate(int lmax, double gm, double r0) {
    ShGrav sh;
    sh.lmax     = lmax;
    sh.cilm_dim = lmax + 1;
    sh.gm       = gm;
    sh.r0       = r0;
    sh.cilm.assign(2 * sh.cilm_dim * sh.cilm_dim, 0.0);
    return sh;
}

int shGravLoadCoefficients(ShGrav &sh, const char *filename) {
    int exitstatus = 0;
    int lmax_read  = sh.lmax;
    double header[8];
    shtools::SHRead(filename,
                    static_cast<int>(strlen(filename)),
                    sh.cilm.data(),
                    sh.cilm_dim,
                    &lmax_read,
                    nullptr,
                    header,
                    8,
                    nullptr,
                    &exitstatus);
    return (exitstatus == 0) ? 1 : -std::abs(exitstatus);
}

void shGravCalculate(const ShGrav &sh,
                     double x, double y, double z,
                     double result[3])
{
    double r   = std::sqrt(x*x + y*y + z*z);
    double lat = std::asin(z / r) * 180.0 / M_PI;
    double lon = std::atan2(y, x) * 180.0 / M_PI;

    // shtools returns (g_r, g_theta, g_phi) with theta = colatitude,
    double res[3];
    shtools::MakeGravGridPoint(sh.cilm.data(), sh.cilm_dim, sh.lmax,
                               sh.gm, sh.r0, r, lat, lon, res,
                               nullptr, nullptr);

    double lat_rad = lat * M_PI / 180.0;
    double lon_rad = lon * M_PI / 180.0;
    double cl = std::cos(lat_rad), sl = std::sin(lat_rad);
    double cp = std::cos(lon_rad), sp = std::sin(lon_rad);

    // e_r       = ( cl*cp,  cl*sp,  sl)
    // e_theta   = ( sl*cp,  sl*sp, -cl)
    // e_phi     = (-sp,     cp,     0)
    result[0] = res[0]*cl*cp + res[1]*sl*cp - res[2]*sp;
    result[1] = res[0]*cl*sp + res[1]*sl*sp + res[2]*cp;
    result[2] = res[0]*sl    - res[1]*cl;
}

void shGravGetPartialsMatrix(const ShGrav &sh,
                             double x, double y, double z,
                             double M[9])
{
    double r = std::sqrt(x*x + y*y + z*z);
    double h = r * 1e-5;

    for (int j = 0; j < 3; j++) {
        double xyz_p[3] = {x, y, z};
        double xyz_m[3] = {x, y, z};
        xyz_p[j] += h;
        xyz_m[j] -= h;

        double g_p[3], g_m[3];
        shGravCalculate(sh, xyz_p[0], xyz_p[1], xyz_p[2], g_p);
        shGravCalculate(sh, xyz_m[0], xyz_m[1], xyz_m[2], g_m);

        for (int i = 0; i < 3; i++)
            M[i*3 + j] = (g_p[i] - g_m[i]) / (2.0 * h);
    }

    for (int i = 0; i < 3; i++)
        for (int j = i+1; j < 3; j++) {
            double s = 0.5 * (M[i*3+j] + M[j*3+i]);
            M[i*3+j] = M[j*3+i] = s;
        }
}

int shGravGetCoefficientC(const ShGrav &sh, int l, int m, double *c) {
    if (l < 0 || l > sh.lmax || m < 0 || m > l) return -1;
    *c = sh.cilm[shGravIndex(0, l, m, sh.cilm_dim)];
    return 1;
}
int shGravGetCoefficientS(const ShGrav &sh, int l, int m, double *s) {
    if (l < 0 || l > sh.lmax || m < 0 || m > l) return -1;
    *s = sh.cilm[shGravIndex(1, l, m, sh.cilm_dim)];
    return 1;
}
int shGravSetCoefficientC(ShGrav &sh, int l, int m, double c) {
    if (l < 0 || l > sh.lmax || m < 0 || m > l) return -1;
    sh.cilm[shGravIndex(0, l, m, sh.cilm_dim)] = c;
    return 1;
}
int shGravSetCoefficientS(ShGrav &sh, int l, int m, double s) {
    if (l < 0 || l > sh.lmax || m < 0 || m > l) return -1;
    sh.cilm[shGravIndex(1, l, m, sh.cilm_dim)] = s;
    return 1;
}