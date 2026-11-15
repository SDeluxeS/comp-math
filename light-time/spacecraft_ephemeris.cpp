#include "include/spacecraft_ephemeris.hpp"

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

    std::vector<double> times;
    std::vector<double> xs;
    std::vector<double> ys;
    std::vector<double> zs;

    std::string line;

    while (std::getline(input, line)) {
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

        while (row >> token) {
            double value = 0.0;
            if (parseNumber(token, value)) {
                values.push_back(value);
            }
        }

        if (values.size() < 4) {
            continue;
        }

        times.push_back(values[0]);
        xs.push_back(values[1] * positionScale);
        ys.push_back(values[2] * positionScale);
        zs.push_back(values[3] * positionScale);
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