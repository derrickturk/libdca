#ifndef BESTFIT_HPP

#include "exponential.hpp"
#include "hyperbolic.hpp"
#include "hyptoexp.hpp"

#include "convex.hpp"
#include "tuple_tools.hpp"

#include <tuple>
#include <numeric>
#include <algorithm>
#include <functional>
#include <cmath>
#include <limits>

namespace dca {

namespace detail {

template<class Decline, class RateIter, class TimeIter>
inline double sse_against_rate(const Decline& decl,
        RateIter rate_begin, RateIter rate_end, TimeIter time_begin)
{
    return std::inner_product(rate_begin, rate_end, time_begin, 0.0, 
            std::plus<double>(),
            [&](double rate, double time) {
                return std::pow(rate - decl.rate(time), 2);
            });
}

template<class Decline, class VolIter>
inline double sse_against_interval(const Decline& decl,
        VolIter vol_begin, VolIter vol_end,
        double time_initial, double time_step)
{
    struct cumulator {
        const Decline& d;
        double cum;
        double last_cum;
        double t;
        double step;

        cumulator(const Decline& d, double t_init, double t_step)
            : d(d), cum(0.0), last_cum(0.0), t(t_init), step(t_step)
        { }

        double operator()(double sse, double vol)
        {
            t += step;
            double interval = d.cumulative(t) - last_cum;
            last_cum += interval;
            return sse + std::pow(vol - interval, 2);
        }
    };

    return std::accumulate(vol_begin, vol_end, 0.0,
            cumulator(decl, time_initial, time_step));
}

template<class Decline>
struct decline_traits {
};

template<>
struct decline_traits<arps_exponential> {
    static convex::simplex<double, double> initial_simplex() noexcept
    {
        return convex::simplex<double, double> {
            std::make_tuple(1, 0.01),
            std::make_tuple(1e4, 5.0),
            std::make_tuple(5e2, 2.3)
        };
    }
};

template<>
struct decline_traits<arps_hyperbolic> {
    static convex::simplex<double, double, double> initial_simplex() noexcept
    {
        return convex::simplex<double, double, double> {
            std::make_tuple(1, 0.01, 0.1),
            std::make_tuple(1e4, 5.0, 5.0),
            std::make_tuple(5e2, 2.3, 2.0),
            std::make_tuple(50, 1.0, 0.75)
        };
    }
};

template<>
struct decline_traits<arps_hyperbolic_to_exponential> {
    static convex::simplex<double, double, double, double> initial_simplex()
      noexcept
    {
        return convex::simplex<double, double, double, double> {
            std::make_tuple(1, 0.01, 0.1, 0.05),
            std::make_tuple(1e4, 5.0, 5.0, 0.05),
            std::make_tuple(5e2, 2.3, 2.0, 0.15),
            std::make_tuple(1e3, 1.5, 1.5, 0.10),
            std::make_tuple(50, 1.0, 0.75, 0.05)
        };
    }
};

}

template<class Decline, class RateIter, class TimeIter>
inline Decline best_from_rate(
        RateIter rate_begin, RateIter rate_end, TimeIter time_begin)
{
    return tuple::construct<Decline>(
            convex::nelder_mead(
              [=](const typename decltype(
                detail::decline_traits<Decline>::initial_simplex()
                )::value_type& t) {
                  try {
                      return detail::sse_against_rate(
                          tuple::construct<Decline>(t),
                          rate_begin, rate_end, time_begin);
                  } catch (...) {
                      return std::numeric_limits<double>::infinity();
                  }
              }, detail::decline_traits<Decline>::initial_simplex(), 300));
}

template<class Decline, class VolIter>
inline Decline best_from_interval_volume(
        VolIter vol_begin, VolIter vol_end,
        double time_initial, double time_step)
{
    return tuple::construct<Decline>(
            convex::nelder_mead(
              [=](const typename decltype(
                detail::decline_traits<Decline>::initial_simplex()
                )::value_type& t) {
                  try {
                      return detail::sse_against_interval(
                          tuple::construct<Decline>(t),
                          vol_begin, vol_end, time_initial, time_step);
                  } catch (...) {
                      return std::numeric_limits<double>::infinity();
                  }
              }, detail::decline_traits<Decline>::initial_simplex(), 300));
}

}

#define BESTFIT_HPP
#endif
