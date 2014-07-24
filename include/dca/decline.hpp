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
inline double convert_decline(double D, double b = 1.0) noexcept;

template<> inline
double convert_decline<nominal, nominal>(double D, double) noexcept
{
    return D;
}

template<> inline
double convert_decline<nominal, secant_effective>(double D, double b) noexcept
{
    return 1.0 - std::pow(b * D + 1.0, -1.0 / b);
}

template<> inline
double convert_decline<nominal, tangent_effective>(double D, double) noexcept
{
    return 1.0 - std::exp(-D);
}

template<> inline
double convert_decline<tangent_effective, nominal>(double D, double) noexcept
{
    return -std::log(1.0 - D);
}

template<> inline
double convert_decline<secant_effective, nominal>(double D, double b) noexcept
{
    return (std::pow(1.0 - D, -b) - 1.0) / b;
}

template<> inline
double convert_decline<secant_effective, tangent_effective>(double D, double b)
  noexcept
{
    return convert_decline<nominal, tangent_effective>(
            convert_decline<secant_effective, nominal>(D),
            b);
}

template<> inline
double convert_decline<tangent_effective, tangent_effective>(double D, double)
  noexcept
{
    return D;
}

template<> inline
double convert_decline<tangent_effective, secant_effective>(double D, double b)
  noexcept
{
    return convert_decline<nominal, secant_effective>(
            convert_decline<tangent_effective, nominal>(D),
            b);
}

template<> inline
double convert_decline<secant_effective, secant_effective>(double D, double)
  noexcept
{
    return D;
}

template<decline_rate Type> inline double decline(double D, double b = 1.0)
  noexcept
{
    return convert_decline<Type, nominal>(D, b);
}

template<class Decline, class OutIter>
inline OutIter interval_volumes(const Decline& decline, OutIter out,
        double time_begin, double time_step, std::size_t n)
{
    double cumulative = decline.cumulative(time_begin);
    while (n--) {
        time_begin += time_step;
        double next = decline.cumulative(time_begin);
        *out++ = next - cumulative;
        cumulative = next;
    }
    return out;
}

template<class Decline>
inline double eur(const Decline& decline, double economic_limit,
        double max_time = std::numeric_limits<double>::infinity(),
        double* time_to_eur = nullptr) noexcept
{
    double t_eur = std::min(time_to_rate(decline, economic_limit), max_time);
    if (time_to_eur) *time_to_eur = t_eur;
    return decline.cumulative(t_eur);
}

template<class Decline>
inline double time_to_rate(const Decline& decline, double rate) noexcept
{
    return std::get<0>(convex::nelder_mead([&](double t) {
                return (t < 0.0)
                  ? std::numeric_limits<double>::infinity()
                  : std::abs(decline.rate(t) - rate);
            }, convex::simplex<double> {
                { std::make_tuple(0.0), std::make_tuple(100) }
            }, 300));
}

template<class Decline>
inline double time_to_cumulative(const Decline& decline, double cum) noexcept
{
    return std::get<0>(convex::nelder_mead([&](double t) {
                return (t < 0.0)
                  ? std::numeric_limits<double>::infinity()
                  : std::abs(decline.cumulative(t) - cum);
            }, convex::simplex<double> {
                { std::make_tuple(0.0), std::make_tuple(100) }
            }, 300));
}

}

#define DECLINE_HPP
#endif
