#ifndef EXPONENTIAL_HPP

#include "decline.hpp"
#include <stdexcept>
#include <cmath>

namespace dca {

class arps_exponential {
    public:
        arps_exponential(double qi, double D);

        const double& qi() const noexcept;
        const double& D() const noexcept;

        double rate(double time) const noexcept;
        double cumulative(double time) const noexcept;

    private:
        double qi_;
        double D_;
};

namespace {
    const double eps = 1e-5;
}

inline arps_exponential::arps_exponential(double qi, double D)
    : qi_(qi), D_(D)
{
    if (qi_ < 0.0)
        throw std::out_of_range("qi must be non-negative.");
    if (D_ < 0.0)
        throw std::out_of_range("D must be non-negative.");
}

inline const double& arps_exponential::qi() const noexcept
{
    return qi_;
}

inline const double& arps_exponential::D() const noexcept
{
    return D_;
}

inline double arps_exponential::rate(double time) const noexcept
{
    if (time < 0.0) return 0.0;
    return qi_ * std::exp(-D_ * time);
}

inline double arps_exponential::cumulative(double time) const noexcept
{
    if (D_ < eps)
        return qi_ * time;
    return qi_ / D_ * (1.0 - std::exp(-D_ * time));
}

}

#define EXPONENTIAL_HPP
#endif
