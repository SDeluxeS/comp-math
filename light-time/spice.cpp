#include "spice.hpp"

#include <SpiceUsr.h>

#include <cmath>
#include <stdexcept>
#include <string>

namespace lt {

void Spice::check(const std::string& operation) {
    if (!failed_c()) {
        return;
    }

    SpiceChar message[1841]{};
    getmsg_c("LONG",
             static_cast<SpiceInt>(sizeof(message)),
             message);
    reset_c();

    throw std::runtime_error(operation + ": " + message);
}

Spice::Spice(const std::string& metaKernel) {
    SpiceChar action[] = "RETURN";

    erract_c("SET", 0, action);
    check("erract_c");

    furnsh_c(metaKernel.c_str());
    check("furnsh_c(" + metaKernel + ")");
}

std::array<double, 6> Spice::state(int targetNaif, double etTdbSeconds) const {
    SpiceDouble stateArray[6]{};
    SpiceDouble lightTime{};

    spkez_c(targetNaif,
            etTdbSeconds,
            "J2000",
            "NONE",
            0,
            stateArray,
            &lightTime);
    check("spkez_c");

    return {
        stateArray[0], stateArray[1], stateArray[2],
        stateArray[3], stateArray[4], stateArray[5]
    };
}

std::array<double, 3> Spice::position(int targetNaif, double etTdbSeconds) const {
    std::array<double, 6> _state = state(targetNaif, etTdbSeconds);
    return {_state[0], _state[1], _state[2]};
}

double Spice::ttMinusTdb(double etTdbSeconds) const {
    SpiceDouble stateArray[6]{};
    SpiceDouble lightTime{};

    spkez_c(1000000001,
            etTdbSeconds,
            "J2000",
            "NONE",
            1000000000,
            stateArray,
            &lightTime);
    check("SPICE conversion TDB -> TT");

    return stateArray[0];
}

double Spice::tdbMinusTtFromTt(double etTtSeconds) const {
    SpiceDouble stateArray[6]{};
    SpiceDouble lightTime{};

    spkez_c(1000000001,
            etTtSeconds,
            "J2000",
            "NONE",
            1000000000,
            stateArray,
            &lightTime);
    check("SPICE conversion TT -> TDB");

    return stateArray[0];
}

double Spice::gravitationalParameter(
    int bodyNaif,
    double fallbackKm3PerSec2) const {

    SpiceInt dimension = 0;
    SpiceDouble values[1]{};

    bodvcd_c(bodyNaif, "GM", 1, &dimension, values);

    if (failed_c()) {
        reset_c();
        return fallbackKm3PerSec2;
    }

    if (dimension >= 1 &&
        std::isfinite(values[0]) &&
        values[0] > 0.0) {
        return values[0];
    }

    return fallbackKm3PerSec2;
}

}