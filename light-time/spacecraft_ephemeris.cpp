#include "include/light-time/spacecraft_ephemeris.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace lt {

namespace {

bool parseNumber(std::string token, double& value) {
    std::replace(token.begin(), token.end(), 'D', 'E');
    std::replace(token.begin(), token.end(), 'd', 'E');

    char* end = nullptr;
    value = std::strtod(token.c_str(), &end);

    return end != token.c_str() &&
           *end == '\0' &&
           std::isfinite(value);
}

}

SpacecraftEphemeris::SpacecraftEphemeris(
    const std::string& fileName,
    int centerNaif,
    double positionScale)
    : centerNaif_(centerNaif) {

    if (!std::isfinite(positionScale)) {
        throw std::invalid_argument(
            "Invalid spacecraft position scale");
    }

    std::ifstream input(fileName);
    if (!input) {
        throw std::runtime_error(
            "Cannot open spacecraft ephemeris: " + fileName);
    }

    // orbit_diff.txt stores time as hours since the RK integration epoch.
    // The first part writes this epoch into the header as
    //   # Epoch_TDB_JD=2455646.830000000000
    double epochTdbJd = 0.0;
    bool hasEpoch = false;
    bool formatDetected = false;
    bool isOrbitDiff = false;

    std::vector<double> times;
    std::vector<double> xs;
    std::vector<double> ys;
    std::vector<double> zs;

    std::string line;
    int lineNumber = 0;

    while (std::getline(input, line)) {
        ++lineNumber;

        const std::string epochKey = "Epoch_TDB_JD=";
        const int epochPos = line.find(epochKey);
        if (epochPos != std::string::npos) {
            const std::string epochText =
                line.substr(epochPos + epochKey.size());
            try {
                std::size_t consumed = 0;
                epochTdbJd = std::stod(epochText, &consumed);
                hasEpoch = consumed > 0 && std::isfinite(epochTdbJd);
            } catch (const std::exception&) {
                throw std::runtime_error(
                    "Invalid Epoch_TDB_JD header in " + fileName);
            }
            continue;
        }

        // Strip comments and accept both whitespace- and comma-separated rows.
        const auto comment = line.find('#');
        if (comment != std::string::npos) {
            line.erase(comment);
        }

        for (char& c : line) {
            if (c == ',') {
                c = ' ';
            }
        }

        std::istringstream row(line);
        std::string token;
        std::vector<double> values;
        bool invalidToken = false;

        while (row >> token) {
            double value = 0.0;
            if (!parseNumber(token, value)) {
                invalidToken = true;
                break;
            }
            values.push_back(value);
        }

        if (values.empty()) {
            continue;
        }
        if (invalidToken) {
            throw std::runtime_error(
                "Invalid numeric value in " + fileName +
                " at line " + std::to_string(lineNumber));
        }

        if (!formatDetected) {
            // orbit_diff rows have 15 columns. A normal ephemeris has 4.
            isOrbitDiff = values.size() >= 15;
            formatDetected = true;
        }

        if (isOrbitDiff) {
            if (values.size() < 15) {
                throw std::runtime_error(
                    "orbit_diff row has fewer than 15 columns at line " +
                    std::to_string(lineNumber));
            }
            if (!hasEpoch) {
                throw std::runtime_error(
                    "orbit_diff.txt has no '# Epoch_TDB_JD=...' header. "
                    "Rebuild and rerun the first part to regenerate it.");
            }

            // Column 1: hours from the integration epoch.
            // Columns 4-6: RK model position relative to Mercury, in km.
            const double tdbSeconds =
                (epochTdbJd - 2451545.0) * 86400.0 + values[0] * 3600.0;
            times.push_back(tdbSeconds);
            xs.push_back(values[3] * positionScale);
            ys.push_back(values[4] * positionScale);
            zs.push_back(values[5] * positionScale);
        } else {
            if (values.size() < 4) {
                throw std::runtime_error(
                    "Ephemeris row has fewer than 4 columns at line " +
                    std::to_string(lineNumber));
            }
            // Native ephemeris format:
            // TDB_seconds_from_J2000 X_km Y_km Z_km
            times.push_back(values[0]);
            xs.push_back(values[1] * positionScale);
            ys.push_back(values[2] * positionScale);
            zs.push_back(values[3] * positionScale);
        }
    }

    if (times.size() < 2) {
        throw std::runtime_error(
            "Spacecraft ephemeris must contain at least two points");
    }

    x_.reset(times, xs);
    y_.reset(times, ys);
    z_.reset(times, zs);
}

std::array<double, 3> SpacecraftEphemeris::relativePosition(
    double tdbSeconds) const {
    return {
        x_(tdbSeconds),
        y_(tdbSeconds),
        z_(tdbSeconds)
    };
}

std::array<double, 3> SpacecraftEphemeris::barycentricPosition(
    double tdbSeconds,
    const Spice& spice) const {

    const std::array<double, 3> relative = relativePosition(tdbSeconds);

    if (centerNaif_ == 0) {
        return relative;
    }

    std::array<double, 3> barycentricPosition;
    vadd(spice.position(centerNaif_, tdbSeconds).data(),
         relative.data(),
         barycentricPosition.data());

         return barycentricPosition;
}

}