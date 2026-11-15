#include "include/eop.hpp"

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

} // namespace

EopTable::EopTable(const std::string& fileName) {
    std::ifstream input(fileName);
    if (!input) {
        throw std::runtime_error("Cannot open EOP file: " + fileName);
    }

    std::vector<double> mjdTai;
    std::array<std::vector<double>, 5> columns;
    std::string line;

    while (std::getline(input, line)) {
        const auto comment = line.find('#');
        if (comment != std::string::npos) {
            line.erase(comment);
        }

        for (char& c : line) {
            if (c == ',' || c == ';') {
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

        // MJD_TAI, pmx, pmy, TAI-UT1, dx, dy.
        if (values.size() < 6) {
            continue;
        }

        mjdTai.push_back(values[0]);

        for (std::size_t i = 0; i < 5; ++i) {
            columns[i].push_back(values[i + 1]);
        }
    }

    if (mjdTai.size() < 2) {
        throw std::runtime_error(
            "EOP file must contain at least two data rows");
    }

    for (std::size_t i = 1; i < mjdTai.size(); ++i) {
        if (!(mjdTai[i] > mjdTai[i - 1])) {
            throw std::runtime_error(
                "EOP MJD values must be strictly increasing");
        }
    }

    // Сдвигаем узлы сплайна к началу таблицы: вместо больших MJD
    // интерполируем по небольшому числу суток от первой строки.
    mjdOrigin_ = mjdTai.front();

    std::vector<double> x;
    x.reserve(mjdTai.size());

    for (double mjd : mjdTai) {
        x.push_back(mjd - mjdOrigin_);
    }

    for (std::size_t i = 0; i < splines_.size(); ++i) {
        splines_[i].reset(x, columns[i]);
    }
}

EopValue EopTable::at(const JulianDate& tai) const {
    // Вычитание MJD origin происходит в dd_real.
    const double x = tai.daysSinceMjd(mjdOrigin_);

    return {
        splines_[0](x),
        splines_[1](x),
        splines_[2](x),
        splines_[3](x),
        splines_[4](x)
    };
}

} // namespace lt