#ifndef DECLINE_TRAITS_HPP

#include "exponential.hpp"
#include "hyperbolic.hpp"
#include "hyptoexp.hpp"

#include "convex"

#include <tuple>

namespace dca {

template<class Decline>
struct decline_traits {
};

template<>
struct decline_traits<arps_exponential> {
    convex::simplex<double, double> initial_simplex() noexcept
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
    convex::simplex<double, double> initial_simplex() noexcept
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
    convex::simplex<double, double> initial_simplex() noexcept
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

#define DECLINE_TRAITS_HPP
#endif
