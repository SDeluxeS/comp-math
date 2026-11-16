#pragma once
#include <math.h>

static inline double vdot(const double a[3], const double b[3]) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

static inline double vnorm(const double a[3]) {
    return sqrt(vdot(a, a));
}

static inline void vadd(const double a[3], const double b[3], double r[3]) {
    r[0] = a[0] + b[0]; 
    r[1] = a[1] + b[1]; 
    r[2] = a[2] + b[2];
}

static inline void vscale(const double a[3], const double b, double r[3]) {
    r[0] = b * a[0]; 
    r[1] = b * a[1]; 
    r[2] = b * a[2];
}

static inline void vsub(const double a[3], const double b[3], double r[3]) {
    r[0] = a[0] - b[0]; 
    r[1] = a[1] - b[1]; 
    r[2] = a[2] - b[2];
}

static inline void cross3(const double a[3], const double b[3], double r[3]) {
    r[0] = a[1] * b[2] - a[2] * b[1];
    r[1] = a[2] * b[0] - a[0] * b[2];
    r[2] = a[0] * b[1] - a[1] * b[0];
}