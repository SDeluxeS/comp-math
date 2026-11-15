#pragma once

#include "JulianDate.hpp"
#include "spice.hpp"

namespace lt {

class TimeScales {
public:
    explicit TimeScales(const Spice& spice)
        : spice_(spice) {}

    JulianDate tdbToTt(const JulianDate& tdb) const;
    JulianDate ttToTdb(const JulianDate& tt) const;
    JulianDate ttToTai(const JulianDate& tt) const;

private:
    const Spice& spice_;
};

}