#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "SpiceUsr.h"     
#include "dynamics.h"
#include "integrator.h"

#define AU_KM      149597870.7        
#define DAY_SEC    86400.0            
#define V_SCALE    (AU_KM / DAY_SEC)  

int main(void)
{
    furnsh_c("kernels/de442t.bsp");     
    furnsh_c("kernels/msgr_final.bsp"); 

    if (failed_c()) {
        fprintf(stderr, "Ошибка при загрузке SPICE-ядер\n");
        reset_c();
        return 1;
    }

    const double JD0      = 2455646.81948;
    const double J2000_JD = 2451545.0;
    
    double et0_sec = (JD0 - J2000_JD) * DAY_SEC;   

    double state0[6], lt;
    spkezr_c("MESSENGER", et0_sec, "J2000", "NONE", "MERCURY", state0, &lt);
    if (failed_c()) {
        fprintf(stderr, "Не удалось получить начальное состояние КА\n");
        reset_c();
        return 1;
    }

    double state0_AU[6];
    for (int k = 0; k < 3; ++k) {
        state0_AU[k]   = state0[k] / AU_KM;
        state0_AU[k+3] = state0[k+3] / V_SCALE;
    }

    double h_sec = 60.0;                  
    double h_days = h_sec / DAY_SEC;      
    
    double T_sec = 24.0 * 3600.0;         
    int n_steps = (int)(T_sec / h_sec);  

    double *ts = malloc((n_steps + 1) * sizeof(double));
    double *ys = malloc((n_steps + 1) * STATE_DIM * sizeof(double));
    
    if (!ts || !ys) { 
        fprintf(stderr, "Ошибка выделения памяти\n"); 
        return 1; 
    }

    printf("Начинаем интегрирование на %.1f часов (шаг %.0f с)...\n", T_sec / 3600.0, h_sec);
    integrate(rhs, et0_sec / DAY_SEC, state0_AU, STATE_DIM, h_days, n_steps, ts, ys);
    printf("Интегрирование завершено успешно.\n");

    FILE *fout = fopen("orbit_diff.txt", "w");
    if (!fout) {
        fprintf(stderr, "Не удалось создать файл orbit_diff.txt\n");
        free(ts); free(ys);
        return 1;
    }
    
    /* Записываем подробный заголовок колонок */
    fprintf(fout, "# 1:t_hours  2:|dr|_km  3:|dv|_m/s  "
                  "4:mx_km  5:my_km  6:mz_km  7:mvx_km/s  8:mvy_km/s  9:mvz_km/s  "
                  "10:rx_km  11:ry_km  12:rz_km  13:rvx_km/s  14:rvy_km/s  15:rvz_km/s\n");

    double max_dr = 0.0, max_dv = 0.0;

    for (int i = 0; i <= n_steps; ++i) {
        double t_sec = ts[i] * DAY_SEC;
        double *model_rel_AU = &ys[i * STATE_DIM]; 

        double model_rel_km[6];
        for (int k = 0; k < 3; ++k) {
            model_rel_km[k]   = model_rel_AU[k] * AU_KM;
            model_rel_km[k+3] = model_rel_AU[k+3] * V_SCALE;
        }

        double ref_rel[6];
        spkezr_c("MESSENGER", t_sec, "J2000", "NONE", "MERCURY", ref_rel, &lt);
        if (failed_c()) { reset_c(); continue; }

        double dr = sqrt(
            pow(model_rel_km[0] - ref_rel[0], 2) +
            pow(model_rel_km[1] - ref_rel[1], 2) +
            pow(model_rel_km[2] - ref_rel[2], 2));
            
        double dv = sqrt(
            pow(model_rel_km[3] - ref_rel[3], 2) +
            pow(model_rel_km[4] - ref_rel[4], 2) +
            pow(model_rel_km[5] - ref_rel[5], 2));

        if (dr > max_dr) max_dr = dr;
        if (dv > max_dv) max_dv = dv;

        /* Записываем все данные в одну строку */
        fprintf(fout, "%.6f  %.6e  %.6e  "
                      "%.6e  %.6e  %.6e  %.6e  %.6e  %.6e  "
                      "%.6e  %.6e  %.6e  %.6e  %.6e  %.6e\n",
                (t_sec - et0_sec) / 3600.0, dr, dv * 1000.0,
                model_rel_km[0], model_rel_km[1], model_rel_km[2], // Модель x,y,z
                model_rel_km[3], model_rel_km[4], model_rel_km[5], // Модель vx,vy,vz
                ref_rel[0], ref_rel[1], ref_rel[2],                // Референс x,y,z
                ref_rel[3], ref_rel[4], ref_rel[5]);               // Референс vx,vy,vz
    }
    
    fclose(fout);

    printf("\nРЕЗУЛЬТАТЫ СРАВНЕНИЯ (за %.1f ч)\n", T_sec / 3600.0);
    printf("Максимальное отклонение по координатам |dr| = %.6e км\n", max_dr);
    printf("Максимальное отклонение по скоростям   |dv| = %.6e м/с\n", max_dv * 1000.0);

    free(ts); 
    free(ys);
    kclear_c();

    return 0;
}