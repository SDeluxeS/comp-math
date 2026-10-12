#include "dynamics.h"
#include "SpiceUsr.h"
#include <math.h>

#define AU_KM      149597870.7        /* 1 АЕ в км */
#define DAY_SEC    86400.0            /* 1 сутки в секундах */
#define V_SCALE    (AU_KM / DAY_SEC)  /* Масштаб скорости (км/с в АЕ/сутки) */

/* Скорость света в АЕ/сутки */
#define C_LIGHT    173.1446326846663

const char *BODY_NAMES[N_BODIES] = {
    "SUN", "MERCURY", "VENUS", "EARTH", "MARS BARYCENTER",
    "JUPITER BARYCENTER", "SATURN BARYCENTER", "MOON",
    "URANUS BARYCENTER", "NEPTUNE BARYCENTER", "PLUTO BARYCENTER"
};

const double BODY_GM[N_BODIES] = {
    2.9591220828411956e-04, /* Sun */
    4.9125001948893182e-11, /* Mercury */
    7.2434523326441187e-10, /* Venus */
    8.8876924467071022e-10, /* Earth */
    9.5495488297258119e-11, /* Mars */
    2.8253458252257917e-07, /* Jupiter */
    8.4597059933762903e-08, /* Saturn */
    1.0931894624024351e-11, /* Moon */
    1.2920265649682399e-08, /* Uranus */
    1.5243573478851939e-08, /* Neptune */
    2.1750964786022513e-12  /* Pluto */
};

static inline double vdot(const double a[3], const double b[3]) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

static inline double vnorm(const double a[3]) {
    return sqrt(vdot(a, a));
}

static inline void vsub(const double a[3], const double b[3], double r[3]) {
    r[0] = a[0] - b[0]; 
    r[1] = a[1] - b[1]; 
    r[2] = a[2] - b[2];
}

static void get_body_states(double et_sec, double X[N_BODIES][3], double V[N_BODIES][3], double A[N_BODIES][3]) {
    double state[6], lt;
    for (int b = 0; b < N_BODIES; ++b) {
        spkezr_c(BODY_NAMES[b], et_sec, "J2000", "NONE", "SSB", state, &lt);
        for (int k = 0; k < 3; ++k) {
            X[b][k] = state[k] / AU_KM;            
            V[b][k] = state[k+3] / V_SCALE;        
        }
    }
    for (int b = 0; b < N_BODIES; ++b) {
        A[b][0] = A[b][1] = A[b][2] = 0.0;
        for (int c = 0; c < N_BODIES; ++c) {
            if (c == b) continue;
            double d_cb[3];
            vsub(X[c], X[b], d_cb);
            double r_bc = vnorm(d_cb);
            double f = BODY_GM[c] / (r_bc * r_bc * r_bc);
            A[b][0] += f * d_cb[0];
            A[b][1] += f * d_cb[1];
            A[b][2] += f * d_cb[2];
        }
    }
}

static void eih_acceleration(const double x_A[3], const double v_A[3],
                             const double X[N_BODIES][3], const double V[N_BODIES][3], const double A[N_BODIES][3],
                             int self_idx, double acc_A[3]) 
{
    const double c2 = C_LIGHT * C_LIGHT;
    acc_A[0] = acc_A[1] = acc_A[2] = 0.0;
    double v_A_sq = vdot(v_A, v_A);

    double sum_C_A = 0.0;
    for (int c = 0; c < N_BODIES; ++c) {
        if (c == self_idx) continue;
        double d_ac[3];
        vsub(X[c], x_A, d_ac);
        sum_C_A += BODY_GM[c] / vnorm(d_ac);
    }

    for (int b = 0; b < N_BODIES; ++b) {
        if (b == self_idx) continue;

        const double *x_B = X[b];
        const double *v_B = V[b];
        const double *a_B = A[b];
        double Gm_B = BODY_GM[b];

        double x_AB[3], x_BA[3];
        vsub(x_A, x_B, x_AB);
        vsub(x_B, x_A, x_BA);

        double r_AB = vnorm(x_AB);
        double r_AB_sq = r_AB * r_AB;

        double n_AB[3] = { x_AB[0]/r_AB, x_AB[1]/r_AB, x_AB[2]/r_AB };
        double n_BA[3] = { x_BA[0]/r_AB, x_BA[1]/r_AB, x_BA[2]/r_AB };

        double v_B_sq = vdot(v_B, v_B);
        double vA_dot_vB = vdot(v_A, v_B);
        double nAB_dot_vB = vdot(n_AB, v_B);

        double sum_C_B = 0.0;
        for (int c = 0; c < N_BODIES; ++c) {
            if (c == b) continue;
            double d_bc[3];
            vsub(X[c], x_B, d_bc);
            sum_C_B += BODY_GM[c] / vnorm(d_bc);
        }

        double xBA_dot_aB = vdot(x_BA, a_B);

        double bracket =
              v_A_sq
            + 2.0 * v_B_sq
            - 4.0 * vA_dot_vB
            - 1.5 * (nAB_dot_vB * nAB_dot_vB)
            - 4.0 * sum_C_A
            - sum_C_B
            + 0.5 * xBA_dot_aB;

        double factor_scalar = (Gm_B / r_AB_sq) * (1.0 + bracket / c2);
        acc_A[0] += factor_scalar * n_BA[0];
        acc_A[1] += factor_scalar * n_BA[1];
        acc_A[2] += factor_scalar * n_BA[2];

        double v_diff_term[3] = {
            4.0 * v_A[0] - 3.0 * v_B[0],
            4.0 * v_A[1] - 3.0 * v_B[1],
            4.0 * v_A[2] - 3.0 * v_B[2]
        };
        double nAB_dot_vdiff = vdot(n_AB, v_diff_term);
        double factor3 = (1.0 / c2) * (Gm_B / r_AB_sq) * nAB_dot_vdiff;

        acc_A[0] += factor3 * (v_A[0] - v_B[0]);
        acc_A[1] += factor3 * (v_A[1] - v_B[1]);
        acc_A[2] += factor3 * (v_A[2] - v_B[2]);

        double factor4 = 3.5 / c2 * (Gm_B / r_AB);
        acc_A[0] += factor4 * a_B[0];
        acc_A[1] += factor4 * a_B[1];
        acc_A[2] += factor4 * a_B[2];
    }
}

void rhs(double et_days, const double y[STATE_DIM], double dydt[STATE_DIM]) {
    double X[N_BODIES][3], V[N_BODIES][3], A[N_BODIES][3];
    
    double et_sec = et_days * DAY_SEC;
    get_body_states(et_sec, X, V, A);

    double r_rel[3] = { y[0], y[1], y[2] };
    double v_rel[3] = { y[3], y[4], y[5] };

    double r_msg[3], v_msg[3];
    for (int k = 0; k < 3; k++) {
        r_msg[k] = X[1][k] + r_rel[k]; /* 1 = MERCURY */
        v_msg[k] = V[1][k] + v_rel[k];
    }

    double a_msg[3];
    eih_acceleration(r_msg, v_msg, X, V, A, -1, a_msg);

    double a_merc[3];
    eih_acceleration(X[1], V[1], X, V, A, 1, a_merc);

    dydt[0] = v_rel[0];  
    dydt[1] = v_rel[1];  
    dydt[2] = v_rel[2];
    dydt[3] = a_msg[0] - a_merc[0];
    dydt[4] = a_msg[1] - a_merc[1];
    dydt[5] = a_msg[2] - a_merc[2];
}