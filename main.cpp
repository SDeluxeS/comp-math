#include "pipeline/orbit_runner.hpp"
#include "pipeline/lt_runner.hpp"

#include <iostream>
#include <string>
#include <stdexcept>

int main(int argc, char** argv) {
    if (argc < 7) {
        std::cerr
            << "Usage:\n"
            << " messenger_pipeline <meta.tm> <eop2.txt> "
            << "<stations.txt> <tx_station> <rx_station> "
            << "<t3_TDB_seconds_from_J2000> [output.txt] [lt_reference.txt]\n";

        return 2;
    }

    try {
        // Параметры для Light-Time
        const std::string metaKernel = argv[1];
        const std::string eopFile = argv[2];
        const std::string stationFile = argv[3];
        const std::string txStationId = argv[4];
        const std::string rxStationId = argv[5];

        const double t3TdbSeconds = std::stod(argv[6]);

        // Файлы промежуточных и итоговых результатов
        const std::string spacecraftFile = "results/orbit_diff.txt";

        // Итоговый файл необязателен
        const std::string outputFile =
            argc > 7 ? argv[7] : "";

        // Референсный файл необязателен
        const std::string referenceFile =
            argc > 8 ? argv[8] : "";

        // Координаты аппарата относительно центра Меркурия
        const int spacecraftCenterNaif = 199;
        const double spacecraftPositionScale = 1.0;

        // Этап 1. Интегрирование орбиты
        std::cout << "\n========== ORBIT MODELING ==========\n";

        const int orbitStatus =
            pipeline::runOrbit(spacecraftFile);

        if (orbitStatus != 0) {
            std::cerr << "Ошибка при интегрировании орбиты.\n";
            return orbitStatus;
        }

        std::cout << "Орбитальное моделирование завершено.\n";

        // Этап 2. Решение задачи Light-Time
        std::cout << "\n========== LIGHT-TIME ==========\n";

        const int lightTimeStatus =
            pipeline::runLightTime(
                metaKernel,
                eopFile,
                stationFile,
                spacecraftFile,
                txStationId,
                rxStationId,
                t3TdbSeconds,
                spacecraftCenterNaif,
                spacecraftPositionScale,
                outputFile
            );

        if (lightTimeStatus != 0) {
            std::cerr << "Ошибка при решении задачи Light-Time.\n";
            return lightTimeStatus;
        }

        if (!referenceFile.empty()) {
            std::cout << "\n========== COMPARISON WITH REFERENCE ==========\n";

            const int compareStatus =
                compareWithReference(
                    metaKernel,
                    eopFile,
                    stationFile,
                    spacecraftFile,
                    referenceFile,
                    spacecraftCenterNaif,
                    spacecraftPositionScale
                );

            return compareStatus;
        }

        std::cout << "\nВсе вычисления успешно завершены.\n";

        return 0;
    }
    catch (const std::exception& error) {
        std::cerr << "ERROR: " << error.what() << '\n';
        return 1;
    }
}