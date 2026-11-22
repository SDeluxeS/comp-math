#pragma once

#include "spacecraft_ephemeris.hpp"
#include "stations.hpp"

#include <vector>

namespace lt {

struct LightTimeLeg {
    double emissionTdbSec{};
    double receptionTdbSec{};
    double geometricSec{};
    double shapiroSec{};
    double totalSec{};
    double residualSec{};
    int iterations{};
};

struct RoundTripResult {
    LightTimeLeg up;
    LightTimeLeg down;
};

class LightTimeSolver {
public:
    LightTimeSolver(const Spice& spice,
                    const SpacecraftEphemeris& spacecraft,
                    const StationTransformer& stationTransformer,
                    double cKmPerSec = 299792.458,
                    double gamma = 1.0,
                    double toleranceSec = 1.0e-7);

    LightTimeLeg solveDownlink(double t3TdbSec,
                               const Station& receiver) const;

    LightTimeLeg solveUplink(double t2TdbSec,
                             const Station& transmitter) const;

    RoundTripResult solveRoundTrip(
        double t3TdbSec,
        const Station& transmitter,
        const Station& receiver) const;

private:
    struct Body {
        int naif{};
        double gmKm3PerSec2{};
    };

    double shapiroDelay(double tA,
                        const std::array<double, 3>& rA,
                        double tB,
                        const std::array<double, 3>& rB) const;

    double legEquation(double tA,
                       const std::array<double, 3>& rA,
                       double tB,
                       const std::array<double, 3>& rB,
                       double* geometric = nullptr,
                       double* shapiro = nullptr) const;

    const Spice& spice_;
    const SpacecraftEphemeris& spacecraft_;
    const StationTransformer& stationTransformer_;

    double cKmPerSec_{};
    double gamma_{};
    double toleranceSec_{};

    Body sun_;
    std::vector<Body> otherBodies_;
};

}