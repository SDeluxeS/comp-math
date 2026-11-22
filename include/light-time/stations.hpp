#pragma once

#include "astro_types.hpp"
#include "eop.hpp"
#include "spice.hpp"
#include "time_scales.hpp"

#include <string>
#include <array>
#include <unordered_map>

namespace lt {

struct Station {
    std::string name;

    // ITRF, км, на эпоху JD 2453371.5.
    std::array<double, 3> referencePositionKm{};

    // Скорость станции в ITRF, км/сутки.
    std::array<double, 3> velocityKmPerDay{};

    std::array<double, 3> positionItrf(double tdbSeconds) const;
};

class StationCatalog {
public:
    // Формат DSN-файла:
    // ID diameter type X(m) Y(m) Z(m) VX(m/year) VY(m/year) VZ(m/year).
    static StationCatalog load(const std::string& fileName);

    // Допустимы, например: "14", "DSS-14", "12*" и "DSS-12*".
    const Station& at(const std::string& id) const;

private:
    static std::string normalizeId(const std::string& id);

    std::unordered_map<std::string, Station> stations_;
};

class StationTransformer {
public:
    StationTransformer(const Spice& spice,
                       const TimeScales& timeScales,
                       const EopTable& eop)
        : spice_(spice),
          timeScales_(timeScales),
          eop_(eop) {}

    // Положение станции относительно SSB в J2000, км.
    std::array<double, 3> barycentric(const Station& station,
                       double tdbSeconds) const;

private:
    const Spice& spice_;
    const TimeScales& timeScales_;
    const EopTable& eop_;
};

} // namespace lt