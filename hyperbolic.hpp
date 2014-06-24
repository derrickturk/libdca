#ifndef HYPERBOLIC_HPP

#include "decline.hpp"
#include "exponential.hpp"
#include <stdexcept>
#include <cmath>

namespace dca {

class arps_hyperbolic {
    public:
        arps_hyperbolic(double qi, double Di, double b);

        const double& qi() const noexcept;
        const double& Di() const noexcept;
        const double& b() const noexcept;

        double rate(double time) const noexcept;
        double cumulative(double time) const noexcept;
        double D(double time) const noexcept;

    private:
        double qi_;
        double Di_;
        double b_;

        double harmonic_rate(double time) const noexcept;
        double harmonic_cumulative(double time) const noexcept;

        static constexpr double eps_ = 1e-5;
};

inline arps_hyperbolic::arps_hyperbolic(double qi, double Di, double b)
    : qi_(qi), Di_(Di), b_(b)
{
    if (qi_ < 0.0)
        throw std::out_of_range("qi must be non-negative.");
    if (Di_ < 0.0)
        throw std::out_of_range("Di must be non-negative.");
    if (b < 0.0)
        throw std::out_of_range("b must be non-negative.");
    if (b > 5.0)
        throw std::out_of_range("b is implausibly high.");
}

inline const double& arps_hyperbolic::qi() const noexcept
{
    return qi_;
}

inline const double& arps_hyperbolic::Di() const noexcept
{
    return Di_;
}

inline const double& arps_hyperbolic::b() const noexcept
{
    return b_;
}

inline double arps_hyperbolic::rate(double time) const noexcept
{
    if (time < 0.0) return 0.0;
    if (b_ < eps_) return arps_exponential(qi_, Di_).rate(time);
    if (std::abs(1.0 - b_) < eps_) return harmonic_rate(time);

    return qi_ * std::pow(1.0 + b_ * Di_ * time, -1.0 / b_);
}

inline double arps_hyperbolic::cumulative(double time) const noexcept
{
    if (time <= 0.0) return 0.0;
    if (Di_ < eps_) return qi_ * time;
    if (b_ < eps_) return arps_exponential(qi_, Di_).cumulative(time);
    if (std::abs(1.0 - b_) < eps_) return harmonic_cumulative(time);

    return qi_ / ((1.0 - b_) * Di_) *
        (1.0 - std::pow(1.0 + b_ * Di_ * time, 1.0 - (1.0 / b_)));
}

inline double arps_hyperbolic::D(double time) const noexcept
{
    return Di_ / (1.0 + b_ * Di_ * time);
}

double arps_hyperbolic::harmonic_rate(double time) const noexcept
{
    return qi_ / (1.0 + Di_ * time);
}

double arps_hyperbolic::harmonic_cumulative(double time) const noexcept
{
    return qi_ / Di_ * std::log(1.0 + Di_ * time);
}

#ifndef DCA_NO_IOSTREAMS
#include <iostream>
inline std::ostream& operator<<(std::ostream& os, const arps_hyperbolic& d)
{
    return os << "<Arps hyperbolic decline: (qi = " << d.qi() << ", Di = "
        << d.Di() << ", b = " << d.b() << ")>";
}
#endif

}

#define HYPERBOLIC_HPP
#endif
