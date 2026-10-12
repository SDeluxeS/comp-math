#ifndef INTEGRATOR_H
#define INTEGRATOR_H

/* Определяем тип "указатель на функцию", чтобы интегратор мог принимать любую правую часть */
typedef void (*rhs_fn)(double t, const double y[], double dydt[]);

/* Функция выполнения одного шага интегратора */
void rk8_step(rhs_fn f, double t, const double y[], int n,
              double h, double y_next[]);

/* Функция выполнения полного интегрирования (в цикле) на n_steps шагов */
void integrate(rhs_fn f, double t0, const double y0[], int n,
               double h, int n_steps,
               double *ts, double *ys);

#endif