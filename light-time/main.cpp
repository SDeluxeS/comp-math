#include "include/eop.hpp"
#include "include/light_time.hpp"
#include "include/spacecraft_ephemeris.hpp"
#include "include/spice.hpp"
#include "include/stations.hpp"
#include "include/time_scales.hpp"

#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

int main(int argc, char** argv) {
    if (argc < 8) {
        std::cerr
            << "Usage:\n"
            << "  lighttime <meta.tm> <eop2.txt> "
            << "<stations.txt> <spacecraft.txt> "
            << "<tx_station> <rx_station> "
            << "<t3_TDB_seconds_from_J2000> "
            << "[spacecraft_center_NAIF=1] "
            << "[spacecraft_position_scale=1]\n";
        return 2;
    }

    try {
        const std::string metaKernel = argv[1];
        const std::string eopFile = argv[2];
        const std::string stationFile = argv[3];
        const std::string spacecraftFile = argv[4];
        const std::string txStationId = argv[5];
        const std::string rxStationId = argv[6];
        const double t3TdbSeconds = std::stod(argv[7]);

        // 1 — барицентр системы Меркурия; 199 — центр Меркурия;
        // 0 — координаты аппарата уже относительно SSB.
        const int spacecraftCenterNaif =
            argc > 8 ? std::stoi(argv[8]) : 1;

        const double spacecraftPositionScale =
            argc > 9 ? std::stod(argv[9]) : 1.0;

        lt::Spice spice(metaKernel);
        lt::TimeScales timeScales(spice);
        lt::EopTable eop(eopFile);

        const lt::StationCatalog stations =
            lt::StationCatalog::load(stationFile);

        const lt::SpacecraftEphemeris spacecraft(
            spacecraftFile,
            spacecraftCenterNaif,
            spacecraftPositionScale);

        const lt::StationTransformer stationTransformer(
            spice, timeScales, eop);

        const lt::LightTimeSolver solver(
            spice, spacecraft, stationTransformer);

        const lt::RoundTripResult result =
            solver.solveRoundTrip(
                t3TdbSeconds,
                stations.at(txStationId),
                stations.at(rxStationId));

        std::cout << std::setprecision(15);
        std::cout << "Times are TDB seconds from J2000\n\n";

        std::cout << "t1 (transmit)   = "
                  << result.up.emissionTdbSec << '\n';
        std::cout << "t2 (spacecraft) = "
                  << result.down.emissionTdbSec << '\n';
        std::cout << "t3 (receive)    = "
                  << t3TdbSeconds << "\n\n";

        std::cout << "Uplink geometric = "
                  << result.up.geometricSec << " s\n";
        std::cout << "Uplink Shapiro   = "
                  << result.up.shapiroSec << " s\n";
        std::cout << "Uplink total     = "
                  << result.up.totalSec << " s\n";
        std::cout << "Uplink residual  = "
                  << result.up.residualSec << " s"
                  << ", iterations = "
                  << result.up.iterations << "\n\n";

        std::cout << "Downlink geometric = "
                  << result.down.geometricSec << " s\n";
        std::cout << "Downlink Shapiro   = "
                  << result.down.shapiroSec << " s\n";
        std::cout << "Downlink total     = "
                  << result.down.totalSec << " s\n";
        std::cout << "Downlink residual  = "
                  << result.down.residualSec << " s"
                  << ", iterations = "
                  << result.down.iterations << "\n\n";

        std::cout << "Round-trip total = "
                  << result.up.totalSec + result.down.totalSec
                  << " s\n";

        return 0;
    } catch (const std::exception& error) {
        std::cerr << "ERROR: " << error.what() << '\n';
        return 1;
    }
}
