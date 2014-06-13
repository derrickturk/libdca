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
            double interval = d.cumulative(t) - last_cum;
            last_cum += interval;
            t += step;

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

    template<class RateIter, class TimeIter>
    static arps_exponential best_from_rate(
            RateIter rate_begin, RateIter rate_end, TimeIter time_begin)
    {
        return tuple::construct<arps_exponential>(
                convex::nelder_mead([=](double qi, double D) {
                    return sse_against_rate(arps_exponential(qi, D),
                        rate_begin, rate_end, time_begin);
                    }, initial_simplex(), 300));
    }

    template<class VolIter>
    static arps_exponential best_from_interval_volume(
            VolIter vol_begin, VolIter vol_end,
            double time_initial, double time_step)
    {
        return tuple::construct<arps_exponential>(
                convex::nelder_mead([=](double qi, double D) {
                    return sse_against_interval(arps_exponential(qi, D),
                        vol_begin, vol_end, time_initial, time_step);
                    }, initial_simplex(), 300));
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

    template<class RateIter, class TimeIter>
    static arps_hyperbolic best_from_rate(
            RateIter rate_begin, RateIter rate_end, TimeIter time_begin)
    {
        return tuple::construct<arps_hyperbolic>(
                convex::nelder_mead([=](double qi, double Di, double b) {
                    return sse_against_rate(arps_hyperbolic(qi, Di, b),
                        rate_begin, rate_end, time_begin);
                    }, initial_simplex(), 300));
    }

    template<class VolIter>
    static arps_hyperbolic best_from_interval_volume(
            VolIter vol_begin, VolIter vol_end,
            double time_initial, double time_step)
    {
        return tuple::construct<arps_hyperbolic>(
                convex::nelder_mead([=](double qi, double Di, double b) {
                    return sse_against_interval(arps_hyperbolic(qi, Di, b),
                        vol_begin, vol_end, time_initial, time_step);
                    }, initial_simplex(), 300));
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

    template<class RateIter, class TimeIter>
    static arps_hyperbolic_to_exponential best_from_rate(
            RateIter rate_begin, RateIter rate_end, TimeIter time_begin)
    {
        return tuple::construct<arps_hyperbolic_to_exponential>(
                convex::nelder_mead([=](double qi, double Di,
                        double b, double Df) {
                    return sse_against_rate(
                        arps_hyperbolic_to_exponential(qi, Di, b, Df),
                        rate_begin, rate_end, time_begin);
                    }, initial_simplex(), 300));
    }

    template<class VolIter>
    static arps_hyperbolic_to_exponential best_from_interval_volume(
            VolIter vol_begin, VolIter vol_end,
            double time_initial, double time_step)
    {
        return tuple::construct<arps_hyperbolic_to_exponential>(
                convex::nelder_mead([=](double qi, double Di,
                        double b, double Df) {
                    return sse_against_interval(
                        arps_hyperbolic_to_exponential(qi, Di, b, Df),
                        vol_begin, vol_end, time_initial, time_step);
                    }, initial_simplex(), 300));
    }
};

}

template<class Decline, class RateIter, class TimeIter>
inline Decline best_from_rate(RateIter rate_begin, RateIter rate_end,
        TimeIter time_begin)
{
    return detail::decline_traits<Decline>::best_from_rate(
            rate_begin, rate_end, time_begin);
}

template<class Decline, class VolIter>
inline Decline best_from_interval_volume(VolIter vol_begin, VolIter vol_end,
        double time_initial, double time_step)
{
    return detail::decline_traits<Decline>::best_from_interval_volume(
            vol_begin, vol_end, time_initial, time_step);
}

}

#define BESTFIT_HPP
#endif
