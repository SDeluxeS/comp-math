#pragma once

#include <vector>

namespace lt {

class CubicSpline {
public:
    CubicSpline() = default;

    CubicSpline(const std::vector<double>& x,
                const std::vector<double>& y);

    void reset(const std::vector<double>& x,
               const std::vector<double>& y);

    double operator()(double xQuery) const;

    bool empty() const {
        return x_.empty();
    }

private:
    std::vector<double> x_;
    std::vector<double> y_;
    std::vector<double> secondDerivatives_;
};

}