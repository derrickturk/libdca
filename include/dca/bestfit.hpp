#ifndef BESTFIT_HPP
#define BESTFIT_HPP

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
            {
                std::make_tuple(1, 0.01),
                std::make_tuple(1e6, 0.5),
                std::make_tuple(1e3, 10.0)
            }
        };
    }

    template<class RateIter>
    static std::pair<std::tuple<double, double>, std::tuple<double, double>>
    parameter_bounds_guess(RateIter rate_begin, RateIter rate_end)
    {
        double peak_rate = *std::max_element(rate_begin, rate_end);
        return std::make_pair(
            std::make_tuple(peak_rate / 0.5, 0.0),
            std::make_tuple(peak_rate * 2.0, 10.0)
        );
    }
};

template<>
struct decline_traits<arps_hyperbolic> {
    static convex::simplex<double, double, double> initial_simplex() noexcept
    {
        return convex::simplex<double, double, double> {
            {
                std::make_tuple(1, 0.01, 0.0),
                std::make_tuple(1e6, 1.0, 0.0),
                std::make_tuple(1e5, 10.0, 0.0),
                std::make_tuple(1e4, 5.0, 3.0)
            }
        };
    }

    template<class RateIter>
    static std::pair<std::tuple<double, double, double>,
        std::tuple<double, double, double>>
    parameter_bounds_guess(RateIter rate_begin, RateIter rate_end)
    {
        double peak_rate = *std::max_element(rate_begin, rate_end);
        return std::make_pair(
            std::make_tuple(peak_rate / 0.5, 0.0, 0.0),
            std::make_tuple(peak_rate * 2.0, 10.0, 3.0)
        );
    }
};

template<>
struct decline_traits<arps_hyperbolic_to_exponential> {
    static convex::simplex<double, double, double, double> initial_simplex()
      noexcept
    {
        return convex::simplex<double, double, double, double> {
            {
                std::make_tuple(1, 0.01, 0.1, 0.05),
                std::make_tuple(1e4, 5.0, 5.0, 0.05),
                std::make_tuple(5e2, 2.3, 2.0, 0.15),
                std::make_tuple(1e3, 1.5, 1.5, 0.10),
                std::make_tuple(50, 1.0, 0.75, 0.05)
            }
        };
    }

    template<class RateIter>
    static std::pair<std::tuple<double, double, double, double>,
        std::tuple<double, double, double, double>>
    parameter_bounds_guess(RateIter rate_begin, RateIter rate_end)
    {
        double peak_rate = *std::max_element(rate_begin, rate_end);
        return std::make_pair(
            std::make_tuple(peak_rate / 0.5, 0.0, 0.0, 0.0),
            std::make_tuple(peak_rate * 2.0, 10.0, 3.0, 10.0)
        );
    }
};

}

template<class Decline, class RateIter, class TimeIter>
inline Decline best_from_rate(
        RateIter rate_begin, RateIter rate_end, TimeIter time_begin)
{
    return tuple::construct<Decline>(
            convex::nelder_mead(
              [=](const auto &t) {
                  try {
                      return detail::sse_against_rate(
                          tuple::construct<Decline>(t),
                          rate_begin, rate_end, time_begin);
                  } catch (...) {
                      return std::numeric_limits<double>::infinity();
                  }
              },
              convex::inner_simplex(detail::decline_traits<Decline>::
                  parameter_bounds_guess(rate_begin, rate_end)),
              300));
}

template<class Decline, class VolIter>
inline Decline best_from_interval_volume(
        VolIter vol_begin, VolIter vol_end,
        double time_initial, double time_step)
{
    return tuple::construct<Decline>(
            convex::nelder_mead(
              [=](const auto &t) {
                  try {
                      return detail::sse_against_interval(
                          tuple::construct<Decline>(t),
                          vol_begin, vol_end, time_initial, time_step);
                  } catch (...) {
                      return std::numeric_limits<double>::infinity();
                  }
              },
              convex::inner_simplex(detail::decline_traits<Decline>::
                  parameter_bounds_guess(vol_begin, vol_end)),
              300));
}

}

#endif
