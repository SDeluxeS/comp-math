#pragma once

#include <array>
#include <string>

namespace lt {

class Spice {
public:
    explicit Spice(const std::string& metaKernel);

    std::array<double, 6> state(int targetNaif, double etTdbSeconds) const;
    std::array<double, 3> position(int targetNaif, double etTdbSeconds) const;

    // Возвращает TT - TDB в секундах.
    double ttMinusTdb(double etTdbSeconds) const;

    // Вход — время TT в секундах от J2000.
    // Возвращает TDB - TT в секундах.
    double tdbMinusTtFromTt(double etTtSeconds) const;

    double gravitationalParameter(
        int bodyNaif,
        double fallbackKm3PerSec2) const;

private:
    static void check(const std::string& operation);
};

}