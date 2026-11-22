#pragma once

#include "cubic_spline.hpp"
#include "spice.hpp"
#include "astro_types.hpp"
#include <array>
#include <string>

namespace lt {

class SpacecraftEphemeris {
public:
    // Формат файла:
    // TDB_seconds_from_J2000 X_km Y_km Z_km
    //
    // centerNaif = 1: координаты относительно барицентра системы Меркурия.
    // centerNaif = 199: координаты относительно центра Меркурия,
    //                   если этот объект есть в kernels.
    // centerNaif = 0: координаты уже относительно SSB.
    SpacecraftEphemeris(const std::string& fileName,
                        int centerNaif,
                        double positionScale = 1.0);

    std::array<double, 3> relativePosition(double tdbSeconds) const;

    std::array<double, 3> barycentricPosition(double tdbSeconds,
                               const Spice& spice) const;

private:
    int centerNaif_{};

    CubicSpline x_;
    CubicSpline y_;
    CubicSpline z_;
};

}