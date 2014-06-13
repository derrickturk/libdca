#ifndef DECLINE_HPP

#include "convex.hpp"
#include <cmath>
#include <limits>

namespace dca {

enum decline_rate {
    nominal,
    tangent_effective,
    secant_effective
};

template<decline_rate From, decline_rate To>
inline constexpr double convert_decline(double D, double b = 1.0) noexcept;

template<> inline constexpr
double convert_decline<nominal, nominal>(double D, double) noexcept
{
    return D;
}

template<> inline constexpr
double convert_decline<nominal, secant_effective>(double D, double b) noexcept
{
    return 1.0 - std::pow(b * D + 1.0, -1.0 / b);
}

template<> inline constexpr
double convert_decline<nominal, tangent_effective>(double D, double) noexcept
{
    return 1.0 - std::exp(-D);
}

template<> inline constexpr
double convert_decline<tangent_effective, nominal>(double D, double) noexcept
{
    return -std::log(1.0 - D);
}

template<> inline constexpr
double convert_decline<secant_effective, nominal>(double D, double b) noexcept
{
    return (std::pow(1.0 - D, -b) - 1.0) / b;
}

template<> inline constexpr
double convert_decline<secant_effective, tangent_effective>(double D, double b)
  noexcept
{
    return convert_decline<nominal, secant_effective>(
            convert_decline<tangent_effective, nominal>(D),
            b);
}

template<> inline constexpr
double convert_decline<tangent_effective, tangent_effective>(double D, double)
  noexcept
{
    return D;
}

template<> inline constexpr
double convert_decline<tangent_effective, secant_effective>(double D, double b)
  noexcept
{
    return convert_decline<nominal, secant_effective>(
            convert_decline<tangent_effective, nominal>(D),
            b);
}

template<> inline constexpr
double convert_decline<secant_effective, secant_effective>(double D, double)
  noexcept
{
    return D;
}

template<decline_rate Type> inline constexpr double decline(double D) noexcept
{
    return convert_decline<Type, nominal>(D);
}

template<class Decline>
inline double eur(const Decline& decline, double economic_limit,
        double max_time = std::numeric_limits<double>::infinity()) noexcept
{
    return decline.cumulative(std::min(time_to_rate(decline, economic_limit),
                max_time));
}

template<class Decline>
inline double time_to_rate(const Decline& decline, double rate) noexcept
{
    return std::get<0>(convex::nelder_mead([&](const std::tuple<double>& t) {
                return (std::get<0>(t) < 0.0)
                  ? std::numeric_limits<double>::infinity()
                  : std::abs(decline.rate(std::get<0>(t)) - rate);
            }, convex::simplex<double> {
                std::make_tuple(0.0), std::make_tuple(100)
            }, 300));
}

template<class Decline>
inline double time_to_cumulative(const Decline& decline, double cum) noexcept
{
    return std::get<0>(convex::nelder_mead([&](const std::tuple<double>& t) {
                return (std::get<0>(t) < 0.0)
                  ? std::numeric_limits<double>::infinity()
                  : std::abs(decline.cumulative(std::get<0>(t)) - cum);
            }, convex::simplex<double> {
                std::make_tuple(0.0), std::make_tuple(100)
            }, 300));
}

}

#define DECLINE_HPP
#endif
