#pragma once

#include "JulianDate.hpp"
#include "cubic_spline.hpp"

#include <array>
#include <string>
#include <vector>

namespace lt {

struct EopValue {
    double pmxMas{};
    double pmyMas{};
    double taiMinusUt1Ms{};
    double dxMas{};
    double dyMas{};
};

class EopTable {
public:
    explicit EopTable(const std::string& fileName);

    // Аргумент — дата TAI. Интерполяция выполняется по MJD в TAI.
    EopValue at(const JulianDate& tai) const;

private:
    double mjdOrigin_{};
    std::array<CubicSpline, 5> splines_;
};

} // namespace lt