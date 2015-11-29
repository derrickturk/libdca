#ifndef CONVEX_HPP
#define CONVEX_HPP

#include <tuple>
#include <utility>
#include <type_traits>
#include <array>
#include <algorithm>
#include <cstddef>
#include <cmath>

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

template<class Fn, class Tuple>
struct must_apply {
    template<class FnDep = Fn>
    static decltype(std::declval<FnDep>()(std::declval<Tuple>()),
            std::false_type {})
    must_apply_test(void*);

    template<class = void>
    static std::true_type must_apply_test(...);

    using type = typename std::decay<decltype(must_apply_test(nullptr))>::type;

    static const bool value = type::value;
};

}

template<class... Params>
using simplex =
    typename detail::simplex_traits<std::tuple<Params...>>::simplex_type;

template<class Fn, class Simplex,
    class = typename std::enable_if<!detail::must_apply<
      Fn, typename Simplex::value_type>::value>::type>
typename Simplex::value_type nelder_mead(
        Fn f,
        const Simplex& initial_simplex,
        int max_iter,
        double term_eps = std::sqrt(std::numeric_limits<double>::epsilon()),
        int term_iter = 10,
        double ref_factor = 1.0,
        double exp_factor = 2.0,
        double con_factor = 0.5,
        double shr_factor = 0.5);

template<class Fn, class Simplex,
    class = typename std::enable_if<detail::must_apply<
      Fn, typename Simplex::value_type>::value>::type,
    class = void>
typename Simplex::value_type nelder_mead(
        Fn f,
        const Simplex& initial_simplex,
        int max_iter,
        double term_eps = std::sqrt(std::numeric_limits<double>::epsilon()),
        int term_iter = 10,
        double ref_factor = 1.0,
        double exp_factor = 2.0,
        double con_factor = 0.5,
        double shr_factor = 0.5);

template<class Fn, class Simplex, class>
typename Simplex::value_type nelder_mead(
        Fn f,
        const Simplex& initial_simplex,
        int max_iter,
        double term_eps, int term_iter,
        double ref_factor,
        double exp_factor,
        double con_factor,
        double shr_factor)
{
    using std::begin;
    using std::end;

    auto trial_simplex(initial_simplex);
    std::array<typename std::result_of<Fn(typename Simplex::value_type)>::type,
        std::tuple_size<Simplex>::value> result;
    std::transform(begin(trial_simplex), end(trial_simplex), begin(result), f);

    auto extrema = std::minmax_element(begin(result), end(result));
    std::size_t best = static_cast<std::size_t>(std::distance(begin(result),
                extrema.first)),
        worst = static_cast<std::size_t>(std::distance(begin(result),
                    extrema.second));
    auto cent = detail::centroid(trial_simplex, worst);

    for (int i = 0, t = 0; t < term_iter && i < max_iter; ++i) {
        auto reflect = detail::tuple_2_scale_add(
                cent, 1.0 + ref_factor, trial_simplex[worst], -ref_factor);
        auto reflect_res = f(reflect);

        if (reflect_res < result[best]) {
            // reflection was better than the best, try expanding
            auto expand = detail::tuple_2_scale_add(
                    cent, 1.0 - exp_factor, reflect, exp_factor);
            auto expand_res = f(expand);
            if (expand_res < reflect_res) {
                trial_simplex[worst] = expand;
                result[worst] = expand_res;
            } else {
                trial_simplex[worst] = reflect;
                result[worst] = reflect_res;
            }

            best = worst;
            worst = static_cast<std::size_t>(std::distance(
                        begin(result),
                        std::max_element(begin(result), end(result))));
            cent = detail::centroid(trial_simplex, worst);
        } else { // reflection was not better than the best
            bool reflection_better_than_second_worst = false;
            for (size_t i = 0; i < result.size(); ++i)
                if (i != worst && result[i] > reflect_res)
                    reflection_better_than_second_worst = true;

            if (reflection_better_than_second_worst) {
                // accept reflected point
                trial_simplex[worst] = reflect;
                result[worst] = reflect_res;
                worst = static_cast<std::size_t>(std::distance(
                            begin(result),
                            std::max_element(begin(result), end(result))));
                cent = detail::centroid(trial_simplex, worst);
            } else {
                // try contracting

                if (result[worst] > reflect_res) {
                    // better than worst: outside contraction
                    auto contract = detail::tuple_2_scale_add(
                            cent, 1.0 - con_factor,
                            reflect, con_factor);
                    auto contract_res = f(contract);
                    if (contract_res <= reflect_res) {
                        trial_simplex[worst] = contract;
                        result[worst] = contract_res;
                        worst = static_cast<std::size_t>(std::distance(
                                    begin(result),
                                    std::max_element(begin(result), end(result))));
                        cent = detail::centroid(trial_simplex, worst);
                    } else { // shrink everything toward best
                        for (std::size_t i = 0; i < trial_simplex.size(); ++i)
                            if (i != best)
                                trial_simplex[i] = detail::tuple_2_scale_add(
                                        trial_simplex[best], 1.0 - shr_factor,
                                        trial_simplex[i], shr_factor);

                        std::transform(begin(trial_simplex), end(trial_simplex),
                                begin(result), f);
                        extrema = std::minmax_element(begin(result), end(result));
                        best = static_cast<std::size_t>(std::distance(
                                    begin(result), extrema.first));
                        worst = static_cast<std::size_t>(std::distance(
                                    begin(result), extrema.second));
                        cent = detail::centroid(trial_simplex, worst);
                    }
                } else { // as bad as worst: inside contraction
                    auto contract = detail::tuple_2_scale_add(
                            cent, 1.0 - con_factor,
                            trial_simplex[worst], con_factor);
                    auto contract_res = f(contract);
                    if (contract_res < result[worst]) {
                        trial_simplex[worst] = contract;
                        result[worst] = contract_res;
                        worst = static_cast<std::size_t>(std::distance(
                                    begin(result),
                                    std::max_element(begin(result), end(result))));
                        cent = detail::centroid(trial_simplex, worst);
                    } else { // shrink everything toward best
                        for (std::size_t i = 0; i < trial_simplex.size(); ++i)
                            if (i != best)
                                trial_simplex[i] = detail::tuple_2_scale_add(
                                        trial_simplex[best], 1.0 - shr_factor,
                                        trial_simplex[i], shr_factor);

                        std::transform(begin(trial_simplex), end(trial_simplex),
                                begin(result), f);
                        extrema = std::minmax_element(begin(result), end(result));
                        best = static_cast<std::size_t>(std::distance(
                                    begin(result), extrema.first));
                        worst = static_cast<std::size_t>(std::distance(
                                    begin(result), extrema.second));
                        cent = detail::centroid(trial_simplex, worst);
                    }
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

template<class Fn, class Simplex, class, class>
typename Simplex::value_type nelder_mead(
        Fn f,
        const Simplex& initial_simplex,
        int max_iter,
        double term_eps, int term_iter,
        double ref_factor,
        double exp_factor,
        double con_factor,
        double shr_factor)
{
    return nelder_mead(
            [&](const typename Simplex::value_type& t) {
                return tuple::apply(f, t);
            },
            initial_simplex,
            max_iter,
            term_eps, term_iter,
            ref_factor, exp_factor, con_factor, shr_factor);
}

}

#endif
