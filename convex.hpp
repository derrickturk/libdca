#ifndef CONVEX_HPP

#include <tuple>
#include <utility>
#include <type_traits>
#include <array>
#include <algorithm>
#include <cstddef>

#include "tuple_tools.hpp"

namespace convex {

namespace detail {

template<class Tuple>
struct simplex_traits {
    using vertex_type = Tuple;
    static const std::size_t vertex_length =
        std::tuple_size<vertex_type>::value;
    static const std::size_t simplex_length = 1 + vertex_length;
    using simplex_type = std::array<vertex_type, simplex_length>;
};

template<std::size_t N, class Tuple>
struct tuple_add_impl {
    static void impl(Tuple& augend, const Tuple& addend)
    {
        std::get<N - 1>(augend) += std::get<N - 1>(addend);
        tuple_add_impl<N - 1, Tuple>::impl(augend, addend);
    }
};

template<class Tuple>
struct tuple_add_impl<0, Tuple> {
    static void impl(Tuple&, const Tuple&)
    {
    }
};

template<class Tuple>
void tuple_add(Tuple& augend,
        const Tuple& addend)
{
    tuple_add_impl<
        std::tuple_size<Tuple>::value,
        Tuple
    >::impl(augend, addend);
}

template<class D, std::size_t N, class Tuple>
struct tuple_divide_scalar_impl {
    static void impl(Tuple& dividend, D divisor)
    {
        std::get<N - 1>(dividend) /= divisor;
        tuple_divide_scalar_impl<D, N - 1, Tuple>::impl(
                dividend, divisor);
    }
};

template<class D, class Tuple>
struct tuple_divide_scalar_impl<D, 0, Tuple> {
    static void impl(Tuple&, D)
    {
    }
};

template<class D, class Tuple>
void tuple_divide_scalar(Tuple& dividend, D divisor)
{
    tuple_divide_scalar_impl<
        D,
        std::tuple_size<Tuple>::value,
        Tuple
    >::impl(dividend, divisor);
}

template<class S, std::size_t N, class Tuple>
struct tuple_2_scale_add_impl {
    static void impl(Tuple& result,
            const Tuple& left, S scale_left,
            const Tuple& right, S scale_right)
    {
        std::get<N - 1>(result) =
            std::get<N - 1>(left) * scale_left +
            std::get<N - 1>(right) * scale_right;
        tuple_2_scale_add_impl<S, N - 1, Tuple>::impl(
                result, left, scale_left, right, scale_right);
    }
};

template<class S, class Tuple>
struct tuple_2_scale_add_impl<S, 0, Tuple> {
    static void impl(Tuple&, const Tuple&, S, const Tuple&, S)
    {
    }
};

template<class S, class Tuple>
Tuple tuple_2_scale_add(const Tuple& left, S scale_left,
        const Tuple& right, S scale_right)
{
    Tuple result;
    tuple_2_scale_add_impl<
        S,
        std::tuple_size<Tuple>::value,
        Tuple
    >::impl(result, left, scale_left, right, scale_right);
    return result;
}

template<class Simplex>
typename Simplex::value_type centroid(const Simplex& spx,
        std::size_t except_index)
{
    typename Simplex::value_type result; // value-initialized (all 0)
    for (std::size_t i = 0; i < spx.size(); ++i) {
        if (i != except_index)
            detail::tuple_add(result, spx[i]);
    }
    detail::tuple_divide_scalar(result, spx.size() - 1);
    return result;
}

}

template<class... Params>
using simplex =
    typename detail::simplex_traits<std::tuple<Params...>>::simplex_type;

template<class Fn, class Simplex>
typename Simplex::value_type nelder_mead(
        Fn f,
        const Simplex& initial_simplex,
        int max_iter,
        double term_eps = std::sqrt(std::numeric_limits<double>::epsilon()),
        int term_iter = 10,
        double ref_factor = 1.0,
        double exp_factor = 2.0,
        double con_factor = 0.5);

template<class Fn, class Simplex>
typename Simplex::value_type nelder_mead(
        Fn f,
        const Simplex& initial_simplex,
        int max_iter,
        double term_eps, int term_iter,
        double ref_factor, double exp_factor, double con_factor)
{
    using std::begin;
    using std::end;

    auto trial_simplex(initial_simplex);
    std::array<typename tuple::result_of<Fn(typename Simplex::value_type)>::type,
        std::tuple_size<Simplex>::value> result;
    std::transform(begin(trial_simplex), end(trial_simplex), begin(result),
            [&](typename Simplex::value_type& t){ return tuple::apply(f, t); });

    auto extrema = std::minmax_element(begin(result), end(result));
    std::size_t best = std::distance(begin(result), extrema.first),
        worst = std::distance(begin(result), extrema.second);
    auto cent = detail::centroid(trial_simplex, worst);

    for (int i = 0, t = 0; t < term_iter && i < max_iter; ++i) {
        auto reflect = detail::tuple_2_scale_add(
                cent, 1.0 + ref_factor, trial_simplex[worst], -ref_factor);
        auto reflect_res = tuple::apply(f, reflect);

        if (reflect_res < result[best]) {
            // reflection was better than the best, try expanding
            auto expand = detail::tuple_2_scale_add(
                    reflect, 1.0 + exp_factor, cent, -exp_factor);
            auto expand_res = tuple::apply(f, expand);
            if (expand_res < result[best]) {
                trial_simplex[worst] = expand;
                result[worst] = expand_res;
            } else {
                trial_simplex[worst] = reflect;
                result[worst] = reflect_res;
            }

            best = worst;
            worst = std::distance(begin(result),
                    std::max_element(begin(result), end(result)));
            cent = detail::centroid(trial_simplex, worst);
        } else {
            auto worse = std::find_if(begin(result), end(result),
                    [=](const decltype(reflect_res)& r) {
                        return r > reflect_res;
                    });

            if (worse != end(result) &&
                    std::distance(begin(result), worse) != worst) {
                // there's somebody worse than the reflected point who is not
                // the worst, so keep the reflected point
                trial_simplex[worst] = reflect;
                result[worst] = reflect_res;
                worst = std::distance(begin(result),
                        std::max_element(begin(result), end(result)));
                cent = detail::centroid(trial_simplex, worst);
            } else {
                // the reflected point is worse than everybody, except
                // for maybe the worst, so contract

                if (worse != end(result)) {
                    // better than the worst!
                    trial_simplex[worst] = reflect;
                    result[worst] = reflect_res;
                    worst = std::distance(begin(result),
                            std::max_element(begin(result), end(result)));
                    cent = detail::centroid(trial_simplex, worst);
                }

                auto contract = detail::tuple_2_scale_add(
                        trial_simplex[worst], con_factor,
                        cent, 1.0 - con_factor);
                auto contract_res = tuple::apply(f, contract);

                if (contract_res >= result[worst]) {
                    // it got worse! (or no better) --- contract everything
                    for (std::size_t i = 0; i < trial_simplex.size(); ++i)
                        if (i != best)
                            trial_simplex[i] = detail::tuple_2_scale_add(
                                    trial_simplex[i], 0.5,
                                    trial_simplex[best], 0.5);

                    std::transform(begin(trial_simplex), end(trial_simplex),
                            begin(result),
                            [&](typename Simplex::value_type& t){
                              return tuple::apply(f, t);
                            });
                    extrema = std::minmax_element(begin(result), end(result));
                    best = std::distance(begin(result), extrema.first);
                    worst = std::distance(begin(result), extrema.second);
                    cent = detail::centroid(trial_simplex, worst);
                } else {
                    trial_simplex[worst] = contract;
                    result[worst] = contract_res;
                    worst = std::distance(begin(result),
                            std::max_element(begin(result), end(result)));
                    cent = detail::centroid(trial_simplex, worst);
                }
            }
        }

        if (result[worst] - result[best] < term_eps)
            ++t;
        else
            t = 0;
    }

    return trial_simplex[best];
}

}

#define CONVEX_HPP
#endif
