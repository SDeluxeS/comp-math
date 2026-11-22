#include "include/light-time/light_time.hpp"
#include "include/body_constants.h"
#include <array>
#include <cmath>
#include <stdexcept>
#include <string>

namespace lt {

namespace {

// Устойчиво вычисляет |a| + |b| - |b-a|.
// Формула через векторное произведение помогает при почти
// противоположных векторах.
double stableSumMinusDistance(const std::array<double, 3>& a,
                              const std::array<double, 3>& b,
                              double normA,
                              double normB,
                              double distanceAB) {
    const double sum = normA + normB;
    const double direct = sum - distanceAB;

    if (vdot(a.data(), b.data()) < 0.0) {
        std::array<double, 3> c{};
        cross3(a.data(), b.data(), c.data());
        const double crossSquared = vdot(c.data(), c.data());

        const double q = normA * normB - vdot(a.data(), b.data());
        const double denominator = q * (sum + distanceAB);

        if (q > 0.0 && denominator > 0.0) {
            const double stable =
                2.0 * crossSquared / denominator;

            if (std::isfinite(stable) && stable > 0.0) {
                return stable;
            }
        }
    }

    return direct;
}

}

LightTimeSolver::LightTimeSolver(
    const Spice& spice,
    const SpacecraftEphemeris& spacecraft,
    const StationTransformer& stationTransformer,
    double cKmPerSec,
    double gamma,
    double toleranceSec)
    : spice_(spice),
      spacecraft_(spacecraft),
      stationTransformer_(stationTransformer),
      cKmPerSec_(cKmPerSec),
      gamma_(gamma),
      toleranceSec_(toleranceSec) {

    if (!(cKmPerSec_ > 0.0) ||
        !(toleranceSec_ > 0.0) ||
        !std::isfinite(gamma_)) {
        throw std::invalid_argument(
            "Invalid light-time solver parameters");
    }

    constexpr double gmScale = BODY_GM_AU3_DAY2_TO_KM3_SEC2;

    sun_ = {10, BODY_GM_SUN_AU3_DAY2 * gmScale};

    struct BodyConstant {
        int naif;
        double gmAu3Day2;
    };

    const std::array<BodyConstant, 10> bodyConstants{{
        {1,   BODY_GM_MERCURY_AU3_DAY2},  // Mercury system
        {2,   BODY_GM_VENUS_AU3_DAY2},    // Venus system
        {399, BODY_GM_EARTH_AU3_DAY2},    // Earth
        {301, BODY_GM_MOON_AU3_DAY2},     // Moon
        {4,   BODY_GM_MARS_AU3_DAY2},     // Mars system
        {5,   BODY_GM_JUPITER_AU3_DAY2},  // Jupiter system
        {6,   BODY_GM_SATURN_AU3_DAY2},   // Saturn system
        {7,   BODY_GM_URANUS_AU3_DAY2},   // Uranus system
        {8,   BODY_GM_NEPTUNE_AU3_DAY2},  // Neptune system
        {9,   BODY_GM_PLUTO_AU3_DAY2}     // Pluto system
    }};

    otherBodies_.reserve(bodyConstants.size());

    for (const BodyConstant& body : bodyConstants) {
        otherBodies_.push_back({
            body.naif,
            body.gmAu3Day2 * gmScale
        });
    }
}

double LightTimeSolver::shapiroDelay(
    double tA,
    const std::array<double, 3>& rA,
    double tB,
    const std::array<double, 3>& rB) const {

    const auto bodyTerm = [&](const Body& body,
                              bool isSun) -> double {
        const std::array<double, 3> bodyA =
            spice_.position(body.naif, tA);
        const std::array<double, 3> bodyB =
            spice_.position(body.naif, tB);

        std::array<double, 3> a{};
        vsub(rA.data(), bodyA.data(), a.data());
        std::array<double, 3> b{};
        vsub(rB.data(), bodyB.data(), b.data());

        const double r1 = vnorm(a.data());
        const double r2 = vnorm(b.data());
        std::array<double, 3> b_diff_a{};
        vsub(b.data(), a.data(), b_diff_a.data());
        const double r12 = vnorm(b_diff_a.data());

        double denominator =
            stableSumMinusDistance(a, b, r1, r2, r12);

        // Солнечное слагаемое уравнения Мойера (8-55)
        // содержит (1+gamma)*GM/c^2 и в числителе,
        // и в знаменателе логарифма.
        if (isSun) {
            denominator +=
                (1.0 + gamma_) * body.gmKm3PerSec2 /
                (cKmPerSec_ * cKmPerSec_);
        }

        if (!(denominator > 0.0) ||
            !std::isfinite(denominator)) {
            throw std::runtime_error(
                "Invalid Shapiro denominator; ray may intersect a body");
        }

        // log((r1+r2+r12+q)/(r1+r2-r12+q))
        // вычисляется через log1p для устойчивости.
        const double logIncrement = 2.0 * r12 / denominator;

        if (!(logIncrement >= 0.0) ||
            !std::isfinite(logIncrement)) {
            throw std::runtime_error(
                "Invalid Shapiro logarithm argument");
        }

        const double logarithm = std::log1p(logIncrement);

        return (1.0 + gamma_) * body.gmKm3PerSec2 /
               (cKmPerSec_ * cKmPerSec_ * cKmPerSec_) *
               logarithm;
    };

    double delay = bodyTerm(sun_, true);

    for (const Body& body : otherBodies_) {
        delay += bodyTerm(body, false);
    }

    return delay;
}

double LightTimeSolver::legEquation(
    double tA,
    const std::array<double, 3>& rA,
    double tB,
    const std::array<double, 3>& rB,
    double* geometric,
    double* shapiro) const {

    std::array<double, 3> delta{};
    vsub(rB.data(), rA.data(), delta.data());

    const double geometricSec = vnorm(delta.data()) / cKmPerSec_;

    const double shapiroSec =
        shapiroDelay(tA, rA, tB, rB);

    if (geometric != nullptr) {
        *geometric = geometricSec;
    }

    if (shapiro != nullptr) {
        *shapiro = shapiroSec;
    }

    return geometricSec + shapiroSec;
}

LightTimeLeg LightTimeSolver::solveDownlink(
    double t3TdbSec,
    const Station& receiver) const {

    const std::array<double, 3> r3 =
        stationTransformer_.barycentric(receiver, t3TdbSec);

    // Начальное приближение: геометрическое время при t2 = t3.
    const std::array<double, 3> r2AtGuess =
        spacecraft_.barycentricPosition(t3TdbSec, spice_);

    std::array<double,3> initialDelta{};  
    vsub(r3.data(), r2AtGuess.data(), initialDelta.data()); 
    double t2 = t3TdbSec - vnorm(initialDelta.data()) / cKmPerSec_;

    double residual = 0.0;
    constexpr int maxIterations = 50;

    for (int iteration = 1;
         iteration <= maxIterations;
         ++iteration) {
        const std::array<double, 3> r2 =
            spacecraft_.barycentricPosition(t2, spice_);

        const double tau =
            legEquation(t2, r2, t3TdbSec, r3);

        const double nextT2 = t3TdbSec - tau;
        const double step = std::abs(nextT2 - t2);

        t2 = nextT2;
        residual = step;

        if (step <= toleranceSec_) {
            const std::array<double, 3> r2Final =
                spacecraft_.barycentricPosition(t2, spice_);

            double geometric = 0.0;
            double shapiro = 0.0;

            const double tauFinal =
                legEquation(t2, r2Final,
                            t3TdbSec, r3,
                            &geometric, &shapiro);

            residual =
                std::abs((t3TdbSec - t2) - tauFinal);

            if (residual <= 2.0 * toleranceSec_) {
                return {
                    t2,
                    t3TdbSec,
                    geometric,
                    shapiro,
                    tauFinal,
                    residual,
                    iteration
                };
            }
        }
    }

    throw std::runtime_error(
        "Downlink light-time solution did not converge; residual=" +
        std::to_string(residual));
}

LightTimeLeg LightTimeSolver::solveUplink(
    double t2TdbSec,
    const Station& transmitter) const {

    const std::array<double, 3> r2 =
        spacecraft_.barycentricPosition(t2TdbSec, spice_);

    const std::array<double, 3> r1AtGuess =
        stationTransformer_.barycentric(transmitter, t2TdbSec);

    std::array<double, 3> initialDelta{};
    vsub(r2.data(), r1AtGuess.data(), initialDelta.data());
    double t1 = t2TdbSec - vnorm(initialDelta.data()) / cKmPerSec_;

    double residual = 0.0;
    constexpr int maxIterations = 50;

    for (int iteration = 1;
         iteration <= maxIterations;
         ++iteration) {
        const std::array<double, 3> r1 =
            stationTransformer_.barycentric(transmitter, t1);

        const double tau =
            legEquation(t1, r1, t2TdbSec, r2);

        const double nextT1 = t2TdbSec - tau;
        const double step = std::abs(nextT1 - t1);

        t1 = nextT1;
        residual = step;

        if (step <= toleranceSec_) {
            const std::array<double, 3> r1Final =
                stationTransformer_.barycentric(transmitter, t1);

            double geometric = 0.0;
            double shapiro = 0.0;

            const double tauFinal =
                legEquation(t1, r1Final,
                            t2TdbSec, r2,
                            &geometric, &shapiro);

            residual =
                std::abs((t2TdbSec - t1) - tauFinal);

            if (residual <= 2.0 * toleranceSec_) {
                return {
                    t1,
                    t2TdbSec,
                    geometric,
                    shapiro,
                    tauFinal,
                    residual,
                    iteration
                };
            }
        }
    }

    throw std::runtime_error(
        "Uplink light-time solution did not converge; residual=" +
        std::to_string(residual));
}

RoundTripResult LightTimeSolver::solveRoundTrip(
    double t3TdbSec,
    const Station& transmitter,
    const Station& receiver) const {

    RoundTripResult result;

    result.down = solveDownlink(t3TdbSec, receiver);
    result.up = solveUplink(result.down.emissionTdbSec, transmitter);

    return result;
}

} // namespace lt