#pragma once

#include <cmath>
#include "include/astro_types.hpp"
#include <qd/dd_real.h>
#include <stdexcept>
#include <array>

namespace lt {

using SofaDate = std::array<double, 2>;
// Julian Date хранится внутри как QD double-double.
// При передаче в SOFA значение разбивается на две части double.
class JulianDate {
public:
    explicit JulianDate(const dd_real& jd)
        : jd_(jd) {}

    static JulianDate fromJulianDate(double jd) {
        return JulianDate(dd_real(jd));
    }

    static JulianDate fromSofaParts(double part1, double part2) {
        return JulianDate(dd_real(part1) + dd_real(part2));
    }

    static JulianDate fromJ2000Seconds(double seconds) {
        if (!std::isfinite(seconds)) {
            throw std::invalid_argument("Non-finite time");
        }

        return fromJ2000Seconds(dd_real(seconds));
    }

    static JulianDate fromJ2000Seconds(const dd_real& seconds) {
        return JulianDate(
            dd_real(2451545.0) +
            seconds / dd_real(86400.0)
        );
    }

    JulianDate addSeconds(double seconds) const {
        if (!std::isfinite(seconds)) {
            throw std::invalid_argument("Non-finite time offset");
        }

        return addSeconds(dd_real(seconds));
    }

    JulianDate addSeconds(const dd_real& seconds) const {
        return JulianDate(
            jd_ + seconds / dd_real(86400.0)
        );
    }

    double daysSince(const JulianDate& other) const {
        return to_double(jd_ - other.jd_);
    }


    double daysSinceJulianDate(double otherJd) const {
        return to_double(jd_ - dd_real(otherJd));
    }

    // MJD минус заданная эпоха MJD.
    // Вычитание выполняется в QD, чтобы не терять точность
    // на больших значениях MJD.
    double daysSinceMjd(double otherMjd) const {
        const dd_real thisMjd = jd_ - dd_real(2400000.5);
        return to_double(thisMjd - dd_real(otherMjd));
    }

    double mjd() const {
        return to_double(jd_ - dd_real(2400000.5));
    }

    double secondsFromJ2000() const {
        return to_double(
            (jd_ - dd_real(2451545.0)) * dd_real(86400.0)
        );
    }

    SofaDate sofaParts() const {
        // Целую часть выбираем через double. Даже если дата очень близка
        // к целому JD и округление перешло через границу, ниже дробная
        // часть нормализуется обратно.
        double part1 = std::floor(to_double(jd_));
        double part2 = to_double(jd_ - dd_real(part1));

        if (part2 < 0.0) {
            part1 -= 1.0;
            part2 += 1.0;
        } else if (part2 >= 1.0) {
            part1 += 1.0;
            part2 -= 1.0;
        }

        return {part1, part2};
    }

private:
    dd_real jd_;
};

}