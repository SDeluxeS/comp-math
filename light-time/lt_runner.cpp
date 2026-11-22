#include "include/pipeline/lt_runner.hpp"
#include "include/light-time/eop.hpp"
#include "include/light-time/light_time.hpp"
#include "include/light-time/spacecraft_ephemeris.hpp"
#include "include/spice.hpp"
#include "include/light-time/stations.hpp"
#include "include/light-time/time_scales.hpp"

#include <iomanip>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

namespace pipeline {}

int runLightTime(
    const std::string& metaKernel,
    const std::string& eopFile,
    const std::string& stationFile,
    const std::string& spacecraftFile,
    const std::string& txStationId,
    const std::string& rxStationId,
    double t3TdbSeconds,
    int spacecraftCenterNaif,
    double spacecraftPositionScale,
    const std::string& outputFile
) {
    try {
        lt::Spice spice(metaKernel);
        lt::TimeScales timeScales(spice);
        lt::EopTable eop(eopFile);

        const lt::StationCatalog stations =
            lt::StationCatalog::load(stationFile);

        const lt::SpacecraftEphemeris spacecraft(
            spacecraftFile,
            spacecraftCenterNaif,
            spacecraftPositionScale
        );

        const lt::StationTransformer stationTransformer(
            spice, timeScales, eop
        );

        const lt::LightTimeSolver solver(
            spice, spacecraft, stationTransformer
        );

        std::cout << "Выполняется расчёт Light-Time...\n";

        const lt::RoundTripResult result =
            solver.solveRoundTrip(
                t3TdbSeconds,
                stations.at(txStationId),
                stations.at(rxStationId)
        );

        auto printResult = [&](std::ostream& out) {
            out << std::setprecision(15);
            out << "Times are TDB seconds from J2000\n\n";

            out << "t1 (transmit)   = "
                    << result.up.emissionTdbSec << '\n';
            out << "t2 (spacecraft) = "
                    << result.down.emissionTdbSec << '\n';
            out << "t3 (receive)    = "
                    << t3TdbSeconds << "\n\n";

            out << "Uplink geometric = "
                    << result.up.geometricSec << " s\n";
            out << "Uplink Shapiro   = "
                    << result.up.shapiroSec << " s\n";
            out << "Uplink total     = "
                    << result.up.totalSec << " s\n";
            out << "Uplink residual  = "
                    << result.up.residualSec << " s"
                    << ", iterations = "
                    << result.up.iterations << "\n\n";

            out << "Downlink geometric = "
                    << result.down.geometricSec << " s\n";
            out << "Downlink Shapiro   = "
                    << result.down.shapiroSec << " s\n";
            out << "Downlink total     = "
                    << result.down.totalSec << " s\n";
            out << "Downlink residual  = "
                    << result.down.residualSec << " s"
                    << ", iterations = "
                    << result.down.iterations << "\n\n";

            out << "Round-trip total = "
                    << result.up.totalSec + result.down.totalSec
                    << " s\n";
        }            

        printResult(std::cout);

        // Дополнительная запись в файл
        if (!outputFile.empty()) {
            std::ofstream output(outputFile);

            if (!output.is_open()) {
                throw std::runtime_error(
                    "Не удалось создать файл результатов: "
                    + outputFile
                );
            }

            printResult(output);

            if (!output) {
                throw std::runtime_error(
                    "Ошибка записи результатов в файл: "
                    + outputFile
                );
            }

            std::cout
                << "\nРезультаты Light-Time записаны в: "
                << outputFile << '\n';
        }

        return 0;
    } catch (const std::exception& error) {
        std::cerr << "ERROR: " << error.what() << '\n';
        return 1;
    }
}

struct ReferenceRecord {
    double t3Jd{};
    double downlinkSec{};
    double uplinkSec{};
    double t1Jd{};
    int dsn{};
};

struct MaxDifference {
    double differenceSec{};
    int lineNumber{};
    int dsn{};
    double referenceValue{};
    double calculatedValue{};
};

std::vector<ReferenceRecord> readReferenceFile(const std::string& referenceFile) {
    std::ifstream input(referenceFile);

    if (!input.is_open())
        throw std::runtime_error("Не удалось открыть референсный файл: " + referenceFile);

    std::vector<ReferenceRecord> records;

    std::string line;
    std::size_t lineNumber = 0;

    while (std::getline(input, line)) {
        ++lineNumber;

        if (line.find_first_not_of(" \t\r\n") == std::string::npos)
            continue;

        // Пропускаем комментарии.
        const std::size_t firstChar = line.find_first_not_of(" \t\r\n");

        if (line[firstChar] == '#') {
            continue;
        }

        std::istringstream stream(line);

        ReferenceRecord record;
        record.lineNumber = lineNumber;

        if (!(stream
            >> record.t3Jd
            >> record.downlinkSec
            >> record.uplinkSec
            >> record.t1Jd
            >> record.dsn
        )) {
            // Некорректная строка или заголовок.
            continue;
        }

        // Проверяем корректность числовых значений.
        if (!std::isfinite(record.t3Jd) ||
            !std::isfinite(record.downlinkSec) ||
            !std::isfinite(record.uplinkSec) ||
            !std::isfinite(record.t1Jd)) {

            throw std::runtime_error("Некорректные числовые значения в строке " + std::to_string(lineNumber));
        }

        if (record.dsn < 0)
            throw std::runtime_error("Некорректный DSN в строке " + std::to_string(lineNumber));


        records.push_back(record);
    }

    if (input.bad())
        throw std::runtime_error("Ошибка чтения референсного файла: " + referenceFile);

    if (records.empty())
        throw std::runtime_error("Референсный файл не содержит корректных записей: " + referenceFile);

    return records;
}

void updateMaximum(MaxDifference& maxDiff, double referenceValue, double calculatedValue, int lineNumber, int dsn) {
    const double diffSec = std::abs(referenceValue - calculatedValue);

    if (diffSec > maxDiff.differenceSec) {
        maxDiff.differenceSec = diffSec;
        maxDiff.lineNumber = lineNumber;
        maxDiff.dsn = dsn;
        maxDiff.referenceValue = referenceValue;
        maxDiff.calculatedValue = calculatedValue;
    }
}

int compareWithReference(
    const std::string& metaKernel,
    const std::string& eopFile,
    const std::string& stationFile,
    const std::string& spacecraftFile,
    const std::string& referenceFile,
    int spacecraftCenterNaif,
    double spacecraftPositionScale
) {
    try {
        std::cout << "\n========================================\n";
        std::cout << "СРАВНЕНИЕ С РЕФЕРЕНСНЫМИ ДАННЫМИ\n";
        std::cout << "========================================\n\n";

        // Читаем референсные данные.
        const std::vector<ReferenceRecord> records = readReferenceFile(referenceFile);

        std::cout << "Количество записей: " << records.size() << '\n';

        // Инициализируем все необходимые компоненты один раз.
        std::cout << "Инициализация SPICE и компонентов модели...\n";

        lt::Spice spice(metaKernel);
        lt::TimeScales timeScales(spice);
        lt::EopTable eop(eopFile);

        const lt::StationCatalog stations =
            lt::StationCatalog::load(stationFile);

        const lt::SpacecraftEphemeris spacecraft(
            spacecraftFile,
            spacecraftCenterNaif,
            spacecraftPositionScale
        );

        const lt::StationTransformer stationTransformer(
            spice, timeScales, eop
        );

        const lt::LightTimeSolver solver(
            spice, spacecraft, stationTransformer
        );

        // Переменные для хранения максимальных расхождений.
        MaxDifference maxDownlink;
        MaxDifference maxUplink;
        MaxDifference maxT1;

        int processedRecords = 0;

        // Обрабатываем каждую строку референса.
        for (const ReferenceRecord& record : records) {

            // Референсный t3 задан в JD TDB.
            // Переводим его в секунды относительно J2000.
            const lt::JulianDate referenceT3 =
                lt::JulianDate::fromJulianDate(record.t3Jd);

            const double t3TdbSeconds = referenceT3.secondsFromJ2000();

            // Станция берётся из столбца DSN.
            const std::string stationId = std::to_string(record.dsn);

            const lt::RoundTripResult result =
                solver.solveRoundTrip(
                    t3TdbSeconds,
                    stations.at(stationId),
                    stations.at(stationId)
                );

            // Сравниваем время распространения сигнала.
            updateMaximum(
                maxDownlink,
                record.downlinkSec,
                result.down.totalSec,
                record.lineNumber,
                record.dsn
            );

            updateMaximum(
                maxUplink,
                record.uplinkSec,
                result.up.totalSec,
                record.lineNumber,
                record.dsn
            );

            // Сравниваем момент отправки сигнала t1.
            // Сначала преобразуем расчётное t1 в JulianDate,
            // затем вычисляем разность через dd_real.

            const lt::JulianDate calculatedT1 =
                lt::JulianDate::fromJ2000Seconds(result.up.emissionTdbSec);

            const lt::JulianDate referenceT1 =
                lt::JulianDate::fromJulianDate(record.t1Jd);

            const double t1DiffSec =
                std::abs(calculatedT1.daysSince(referenceT1) * 86400.0);

            if (t1DiffSec > maxT1.differenceSec) {
                maxT1.differenceSec = t1DiffSec;
                maxT1.lineNumber = record.lineNumber;
                maxT1.dsn = record.dsn;

                maxT1.referenceValue =
                    referenceT1.secondsFromJ2000();

                maxT1.calculatedValue =
                    calculatedT1.secondsFromJ2000();
            }

            ++processedRecords;

            // Периодически выводим прогресс.
            if (processedRecords % 1000 == 0 ||
                processedRecords == records.size()) {

                std::cout << "\rОбработано записей: "
                          << processedRecords << " / "
                          << records.size()
                          << std::flush;
            }
        }

        std::cout << "\n\n";

        // 7. Выводим итоговые максимальные расхождения.
        std::cout << std::setprecision(15);

        std::cout << "========================================\n";
        std::cout << "РЕЗУЛЬТАТЫ СРАВНЕНИЯ\n";
        std::cout << "========================================\n\n";

        std::cout 
            << "Максимальное расхождение DOWNLINK: "
            << maxDownlink << "\n\n"
            << "Максимальное расхождение UPLINK: "
            << maxUplink << "\n\n"
            << "Максимальное расхождение t1: "
            << maxT1 << "\n\n";

        return 0;
    }
    catch (const std::exception& error) {
        std::cerr << "\nERROR при сравнении с референсом: "
                  << error.what() << '\n';

        return 1;
    }
}