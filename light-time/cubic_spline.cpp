#include "cubic_spline.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace lt {

CubicSpline::CubicSpline(const std::vector<double>& x,
                         const std::vector<double>& y) {
    reset(x, y);
}

void CubicSpline::reset(const std::vector<double>& x,
                        const std::vector<double>& y) {
    if (x.size() != y.size() || x.size() < 2) {
        throw std::invalid_argument(
            "CubicSpline requires at least two points");
    }

    for (std::size_t i = 0; i < x.size(); ++i) {
        if (!std::isfinite(x[i]) || !std::isfinite(y[i])) {
            throw std::invalid_argument(
                "CubicSpline input contains non-finite values");
        }

        if (i > 0 && !(x[i] > x[i - 1])) {
            throw std::invalid_argument(
                "CubicSpline x values must be strictly increasing");
        }
    }

    x_ = x;
    y_ = y;
    secondDerivatives_.assign(x.size(), 0.0);

    const std::size_t n = x.size();

    // Естественный сплайн: вторая производная на концах равна нулю.
    if (n == 2) {
        return;
    }

    const std::size_t m = n - 2;

    std::vector<double> lower(m);
    std::vector<double> diagonal(m);
    std::vector<double> upper(m);
    std::vector<double> rhs(m);

    for (std::size_t j = 0; j < m; ++j) {
        const std::size_t i = j + 1;

        const double h0 = x_[i] - x_[i - 1];
        const double h1 = x_[i + 1] - x_[i];

        lower[j] = h0;
        diagonal[j] = 2.0 * (h0 + h1);
        upper[j] = h1;

        rhs[j] = 6.0 * (
            (y_[i + 1] - y_[i]) / h1 -
            (y_[i] - y_[i - 1]) / h0
        );
    }

    // Прямой ход метода прогонки.
    for (std::size_t i = 1; i < m; ++i) {
        const double factor = lower[i] / diagonal[i - 1];
        diagonal[i] -= factor * upper[i - 1];
        rhs[i] -= factor * rhs[i - 1];
    }

    // Обратный ход.
    std::vector<double> solution(m);
    solution[m - 1] = rhs[m - 1] / diagonal[m - 1];

    for (std::size_t j = m - 1; j > 0; --j) {
        solution[j - 1] =
            (rhs[j - 1] - upper[j - 1] * solution[j]) /
            diagonal[j - 1];
    }

    for (std::size_t j = 0; j < m; ++j) {
        secondDerivatives_[j + 1] = solution[j];
    }
}

double CubicSpline::operator()(double xQuery) const {
    if (x_.size() < 2) {
        throw std::runtime_error("CubicSpline is not initialized");
    }

    if (!std::isfinite(xQuery) ||
        xQuery < x_.front() ||
        xQuery > x_.back()) {
        throw std::out_of_range(
            "CubicSpline query is outside the data range");
    }

    auto it = std::upper_bound(x_.begin(), x_.end(), xQuery);

    std::size_t i;
    if (it == x_.begin()) {
        i = 0;
    } else if (it == x_.end()) {
        i = x_.size() - 2;
    } else {
        i = static_cast<std::size_t>(it - x_.begin() - 1);
    }

    const double h = x_[i + 1] - x_[i];
    const double a = (x_[i + 1] - xQuery) / h;
    const double b = (xQuery - x_[i]) / h;

    return a * y_[i] + b * y_[i + 1]
         + ((a * a * a - a) * secondDerivatives_[i]
          + (b * b * b - b) * secondDerivatives_[i + 1])
            * h * h / 6.0;
}

} // namespace lt