#include "include/light-time/stations.hpp"
#include "include/astro_types.hpp"
#include <sofa.h>
#include <cctype>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>

namespace lt {

namespace {

constexpr double kMetersToKm = 1.0e-3;
constexpr double kDaysPerYear = 365.25;

constexpr double kMasToRadians =
    3.14159265358979323846 / (180.0 * 3600.0 * 1000.0);

} // namespace

std::array<double, 3> Station::positionItrf(double tdbSeconds) const {
    const JulianDate tdb =
        JulianDate::fromJ2000Seconds(tdbSeconds);

    const JulianDate referenceEpoch =
        JulianDate::fromJulianDate(2453371.5);

    const double elapsedDays = tdb.daysSince(referenceEpoch);
    std::array<double, 3> result{};
    std::array<double, 3> temp{};
    vscale(velocityKmPerDay.data(), elapsedDays, temp.data());
    vadd(referencePositionKm.data(),temp.data(), result.data());
    return result;
}

std::string StationCatalog::normalizeId(const std::string& id) {
    std::string key;
    key.reserve(id.size());

    for (unsigned char c : id) {
        if (!std::isspace(c)) {
            key.push_back(static_cast<char>(std::toupper(c)));
        }
    }

    // DSS-14, DSS14 и DSS 14 становятся "14".
    // Звёздочка у исторических позиций сохраняется.
    if (key.size() >= 3 && key.compare(0, 3, "DSS") == 0) {
        key.erase(0, 3);

        if (!key.empty() && key.front() == '-') {
            key.erase(0, 1);
        }
    }

    if (key.empty()) {
        throw std::invalid_argument("Empty DSN station ID");
    }

    return key;
}

StationCatalog StationCatalog::load(const std::string& fileName) {
    std::ifstream input(fileName);
    if (!input) {
        throw std::runtime_error(
            "Cannot open station file: " + fileName);
    }

    StationCatalog catalog;
    std::string line;

    while (std::getline(input, line)) {
        const auto comment = line.find('#');
        if (comment != std::string::npos) {
            line.erase(comment);
        }

        std::istringstream row(line);

        std::string stationId;
        std::string diameter;
        std::string antennaType;

        double xM = 0.0;
        double yM = 0.0;
        double zM = 0.0;
        double vxMPerYear = 0.0;
        double vyMPerYear = 0.0;
        double vzMPerYear = 0.0;

        // Например:
        // 14 70m A-E -2353621.3988 ... -0.0156 ...
        if (!(row >> stationId >> diameter >> antennaType
                  >> xM >> yM >> zM
                  >> vxMPerYear >> vyMPerYear >> vzMPerYear)) {
            continue;
        }

        if (!std::isfinite(xM) ||
            !std::isfinite(yM) ||
            !std::isfinite(zM) ||
            !std::isfinite(vxMPerYear) ||
            !std::isfinite(vyMPerYear) ||
            !std::isfinite(vzMPerYear)) {
            throw std::runtime_error(
                "Non-finite data for station " + stationId);
        }

        Station station;
        station.name = stationId;

        // Координаты из файла: m -> km.
        station.referencePositionKm = {
            xM * kMetersToKm,
            yM * kMetersToKm,
            zM * kMetersToKm
        };

        // Скорости из файла: m/year -> km/day.
        const double velocityScale =
            kMetersToKm / kDaysPerYear;

        station.velocityKmPerDay = {
            vxMPerYear * velocityScale,
            vyMPerYear * velocityScale,
            vzMPerYear * velocityScale
        };

        const std::string key = normalizeId(stationId);
        const auto inserted = catalog.stations_.emplace(key, station);

        if (!inserted.second) {
            throw std::runtime_error(
                "Duplicate station ID: " + stationId);
        }
    }

    if (catalog.stations_.empty()) {
        throw std::runtime_error(
            "No station records found in " + fileName);
    }

    return catalog;
}

const Station& StationCatalog::at(const std::string& id) const {
    const std::string key = normalizeId(id);
    const auto it = stations_.find(key);

    if (it == stations_.end()) {
        throw std::out_of_range("Unknown DSN station: " + id);
    }

    return it->second;
}

std::array<double, 3> StationTransformer::barycentric(
    const Station& station,
    double tdbSeconds) const {

    const JulianDate tdb =
        JulianDate::fromJ2000Seconds(tdbSeconds);

    const JulianDate tt = timeScales_.tdbToTt(tdb);
    const JulianDate tai = timeScales_.ttToTai(tt);

    const EopValue eop = eop_.at(tai);

    // UT1 = TAI - (TAI - UT1).
    const double taiMinusUt1Seconds =
        eop.taiMinusUt1Ms / 1000.0;
    const JulianDate ut1 =
        tai.addSeconds(-taiMinusUt1Seconds);

    const SofaDate ttParts = tt.sofaParts();
    const SofaDate ut1Parts = ut1.sofaParts();

    double x = 0.0;
    double y = 0.0;

    iauXy06(ttParts[0], ttParts[1], &x, &y);

    // dx, dy из EOP заданы в milliarcseconds.
    x += eop.dxMas * kMasToRadians;
    y += eop.dyMas * kMasToRadians;

    // Полюсное движение pmx, pmy также задано в milliarcseconds.
    const double xp = eop.pmxMas * kMasToRadians;
    const double yp = eop.pmyMas * kMasToRadians;

    double c2t[3][3]{};

    iauC2txy(ttParts[0], ttParts[1],
             ut1Parts[0], ut1Parts[1],
             x, y, xp, yp, c2t);

    // C2T переводит celestial -> terrestrial.
    // Поэтому для ITRF -> celestial берём транспонированную матрицу.
    const std::array<double, 3> rItrf = station.positionItrf(tdbSeconds);
    std::array<double, 3> stationFromEarth{0.0, 0.0, 0.0};

    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            stationFromEarth[i] += c2t[j][i] * rItrf[j];
        }
    }

    const std::array<double, 3> earthFromSsb =
        spice_.position(399, tdbSeconds);

    std::array<double, 3> result{};
    vadd(earthFromSsb.data(), stationFromEarth.data(), result.data());
    return result;
}

}