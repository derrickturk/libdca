#ifndef DECLINE_HPP

#include <cmath>

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

}

#define DECLINE_HPP
#endif
