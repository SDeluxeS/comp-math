#ifndef DYNAMICS_H
#define DYNAMICS_H

/* Размерность фазового вектора y: [x, y, z, vx, vy, vz] */
#define STATE_DIM 6

/* Количество возмущающих тел Солнечной системы */
#define N_BODIES  11

/* Объявление глобальных массивов */
extern const char *BODY_NAMES[N_BODIES];
extern const double BODY_GM[N_BODIES];

/* 
 * Функция правой части дифференциального уравнения
 */
void rhs(double et, const double y[STATE_DIM], double dydt[STATE_DIM]);

#endif