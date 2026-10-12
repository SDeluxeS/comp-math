#include "integrator.h"
#include "rk8_coeffs.h"
#include <string.h>

void rk8_step(rhs_fn f, double t, const double y[], int n,
              double h, double y_next[])
{
    double k[RK8_STAGES][STATE_DIM_MAX]; 
    double yi[STATE_DIM_MAX];            

    
    for (int i = 0; i < RK8_STAGES; ++i) {
        for (int d = 0; d < n; ++d) {
            yi[d] = y[d];
        }
        for (int j = 0; j < i; ++j) {
            for (int d = 0; d < n; ++d) {
                yi[d] += h * RK8_A[i][j] * k[j][d];
            }
        }
        f(t + RK8_C[i] * h, yi, k[i]);
    }

    for (int d = 0; d < n; ++d) {
        y_next[d] = y[d];
    }
    for (int i = 0; i < RK8_STAGES; ++i) {
        for (int d = 0; d < n; ++d) {
            y_next[d] += h * RK8_B[i] * k[i][d];
        }
    }
}

void integrate(rhs_fn f, double t0, const double y0[], int n,
               double h, int n_steps,
               double *ts, double *ys)
{
    ts[0] = t0;
    memcpy(&ys[0], y0, n * sizeof(double));

    for (int i = 0; i < n_steps; ++i) {
        rk8_step(f, ts[i], &ys[i * n], n, h, &ys[(i + 1) * n]);
        ts[i + 1] = ts[i] + h;
    }
}