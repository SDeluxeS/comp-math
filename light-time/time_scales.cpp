#include "include/time_scales.hpp"
#include <include/sofa.h>

#include <cmath>
#include <stdexcept>

namespace lt {

JulianDate TimeScales::tdbToTt(const JulianDate& tdb) const {
    const double etTdb = tdb.secondsFromJ2000();
    const double offsetSeconds = spice_.ttMinusTdb(etTdb);

    if (!std::isfinite(offsetSeconds)) {
        throw std::runtime_error("Invalid TT-TDB offset from SPICE");
    }

    return tdb.addSeconds(offsetSeconds);
}

JulianDate TimeScales::ttToTdb(const JulianDate& tt) const {
    const double etTt = tt.secondsFromJ2000();
    const double offsetSeconds = spice_.tdbMinusTtFromTt(etTt);

    if (!std::isfinite(offsetSeconds)) {
        throw std::runtime_error("Invalid TDB-TT offset from SPICE");
    }

    return tt.addSeconds(offsetSeconds);
}

JulianDate TimeScales::ttToTai(const JulianDate& tt) const {
    const auto ttParts = tt.sofaParts();

    double tai1 = 0.0;
    double tai2 = 0.0;

    const int status = iauTttai(
        ttParts[0], ttParts[1], &tai1, &tai2);

    if (status < 0) {
        throw std::runtime_error("iauTttai failed");
    }

    return JulianDate::fromSofaParts(tai1, tai2);
}

}