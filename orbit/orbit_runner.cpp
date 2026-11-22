#include <iostream>
#include <cstdio>
#include <string>
#include <math.h>
#include "SpiceUsr.h"     
#include "include/orbit/dynamics.h"
#include "include/orbit/integrator.h"
#include "include/orbit/mercury_gravity.h"
#include "include/pipeline/orbit_runner.hpp"

#define AU_KM      149597870.7        
#define DAY_SEC    86400.0            
#define V_SCALE    (AU_KM / DAY_SEC)  

namespace pipeline {
    int runOrbit(const std::string& outputFile) {
        const double JD0      = 2455642.75;
        const double J2000_JD = 2451545.0;

        double h_sec = 60.0;                  
        double h_days = h_sec / DAY_SEC;      
        
        double T_sec = 26.0 * 3600.0;         
        int n_steps = (int)(T_sec / h_sec);  

        double* ts = nullptr;
        double* ys = nullptr;
        FILE* fout = nullptr;

        bool spiceInitialized = false;

        auto cleanup = [&]() {
            if (fout) {
                fclose(fout);
                fout = nullptr;
            }

            free(ts);
            free(ys);

            ts = nullptr;
            ys = nullptr;

            if (spiceInitialized) {
                kclear_c();
                spiceInitialized = false;
            }
        };

        // Загрузка SPICE-ядер
        furnsh_c("kernels/de442t.bsp");     
        furnsh_c("kernels/msgr_final.bsp"); 
        spiceInitialized = true;

        if (failed_c()) {
            std::cerr << "Ошибка при загрузке SPICE-ядер\n";
            reset_c();
            cleanup();
            return 1;
        }

        /* Инициализация поля HgM009 до начала интегрирования. */
        int grav_status = mercury_gravity_init("kernels/HgM009.sha");

        if (grav_status <= 0) {
            std::cerr << "Ошибка загрузки коэффициентов HgM009.sha (код "
                    << grav_status << ")\n";

            cleanup();
            return 1;
        }

        // Определение начального состояние аппарата
        double et0_sec = (JD0 - J2000_JD) * DAY_SEC;   

        double state0[6], lt;
        spkezr_c("MESSENGER", et0_sec, "J2000", "NONE", "MERCURY", state0, &lt);

        if (failed_c()) {
            std::cerr << "Не удалось получить начальное состояние КА\n";
            cleanup();
            return 1;
        }

        // Перевод начального состояния в АЕ и АЕ/сутки
        double state0_AU[6];
        for (int k = 0; k < 3; ++k) {
            state0_AU[k]   = state0[k] / AU_KM;
            state0_AU[k+3] = state0[k+3] / V_SCALE;
        }


        ts = static_cast<double*>(malloc((n_steps + 1) * sizeof(double)));
        ys = static_cast<double*>(malloc((n_steps + 1) * STATE_DIM * sizeof(double)));
        
        if (!ts || !ys) { 
            std::cerr << "Ошибка выделения памяти\n";
            cleanup();
            return 1; 
        }

        printf("Начинаем интегрирование на %.1f часов (шаг %.0f с)...\n", T_sec / 3600.0, h_sec);
        integrate(rhs, et0_sec / DAY_SEC, state0_AU, STATE_DIM, h_days, n_steps, ts, ys);
        printf("Интегрирование завершено успешно.\n");

        fout = fopen(outputFile.c_str(), "w"); 

        if (!fout) {
            std::cerr << "Не удалось создать файл orbit_diff.txt\n";
            cleanup();
            return 1;
        }
        
        /* Записываем эпоху интегрирования и подробный заголовок колонок.
        Light-Time читает orbit_diff.txt напрямую, без промежуточной конвертации. */
        fprintf(fout, "# Epoch_TDB_JD=%.12f\n", JD0);
        fprintf(fout, "# 1:t_hours  2:|dr|_km  3:|dv|_m/s  "
                    "4:mx_km  5:my_km  6:mz_km  7:mvx_km/s  8:mvy_km/s  9:mvz_km/s  "
                    "10:rx_km  11:ry_km  12:rz_km  13:rvx_km/s  14:rvy_km/s  15:rvz_km/s\n");

        double max_dr = 0.0, max_dv = 0.0;

        for (int i = 0; i <= n_steps; ++i) {
            double t_sec = ts[i] * DAY_SEC;
            double* model_rel_AU = &ys[i * STATE_DIM]; 

            double model_rel_km[6];
            for (int k = 0; k < 3; ++k) {
                model_rel_km[k]   = model_rel_AU[k] * AU_KM;
                model_rel_km[k+3] = model_rel_AU[k+3] * V_SCALE;
            }

            double ref_rel[6];
            spkezr_c("MESSENGER", t_sec, "J2000", "NONE", "MERCURY", ref_rel, &lt);

            if (failed_c()) { reset_c(); cleanup(); continue; }

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
            fprintf(fout, "%.12f  %.6e  %.6e  "
                        "%.6e  %.6e  %.6e  %.6e  %.6e  %.6e  "
                        "%.6e  %.6e  %.6e  %.6e  %.6e  %.6e\n",
                    (t_sec - et0_sec) / 3600.0, dr, dv * 1000.0,
                    model_rel_km[0], model_rel_km[1], model_rel_km[2], // Модель x,y,z
                    model_rel_km[3], model_rel_km[4], model_rel_km[5], // Модель vx,vy,vz
                    ref_rel[0], ref_rel[1], ref_rel[2],                // Референс x,y,z
                    ref_rel[3], ref_rel[4], ref_rel[5]);               // Референс vx,vy,vz
        }
        

        printf("\nРЕЗУЛЬТАТЫ СРАВНЕНИЯ (за %.1f ч)\n", T_sec / 3600.0);
        printf("Максимальное отклонение по координатам |dr| = %.6e км\n", max_dr);
        printf("Максимальное отклонение по скоростям   |dv| = %.6e м/с\n", max_dv * 1000.0);

        cleanup();

        return 0;
    }
}