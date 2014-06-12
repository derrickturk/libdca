#ifndef CONVEX_HPP

#include <utility>
#include <type_traits>
#include <array>
#include <algorithm>
#include <cstddef>

namespace convex {

namespace detail {

template<class... Params>
struct simplex_traits {
    using vertex_type = std::tuple<Params...>;
    static const std::size_t vertex_length =
        std::tuple_size<vertex_type>::value;
    static const std::size_t simplex_length = 1 + vertex_length;
    using simplex_type = std::array<vertex_type, simplex_length>;
};

template<std::size_t N, class... Params>
struct tuple_add_impl {
    static void impl(std::tuple<Params...>& augend,
            const std::tuple<Params...>& addend)
    {
        std::get<N - 1>(augend) += std::get<N - 1>(addend);
        typename tuple_add_impl<N - 1, Params...>::impl(
                augend, addend);
    }
};

template<class... Params>
struct tuple_add_impl<0, Params...> {
    static void impl(std::tuple<Params...>& augend,
            const std::tuple<Params...>& addend)
    {
    }
};

template<class... Params>
void tuple_add(std::tuple<Params...>& augend,
        const std::tuple<Params...>& addend)
{
    tuple_add_impl<
        std::tuple_size<std::tuple<Params...>>::value,
        Params...
    >::impl(augend, addend);
}

template<class D, std::size_t N, class... Params>
struct tuple_divide_scalar_impl {
    static void impl(std::tuple<Params...>& dividend, D divisor)
    {
        std::get<N - 1>(dividend) /= divisor;
        typename tuple_divide_scalar_impl<D, N - 1, Params...>::impl(
                dividend, divisor);
    }
};

template<class D, class... Params>
struct tuple_divide_scalar_impl<D, 0, Params...> {
    static void impl(std::tuple<Params...>& dividend, D divisor)
    {
    }
};

template<class D, class... Params>
void tuple_divide_scalar(std::tuple<Params...>& dividend, D divisor)
{
    tuple_divide_scalar_impl<
        D,
        std::tuple_size<std::tuple<Params...>>::value,
        Params...
    >::impl(dividend, divisor);
}

}

template<class... Params>
using simplex = typename detail::simplex_traits<Params...>::simplex_type;

template<class... Params>
std::tuple<Params...> centroid(const simplex<Params...>& spx,
        std::size_t except_index);

template<class Fn, class Simplex>
typename Simplex::value_type nelder_mead(
        Fn f,
        const Simplex& initial_simplex,
        int max_iter,
        double ref_factor = 1.0,
        double exp_factor = 2.0,
        double con_factor = 0.5);

template<class... Params>
std::tuple<Params...> centroid(const simplex<Params...>& spx,
        std::size_t except_index)
{
    std::tuple<Params...> result; // value-initialized (all 0)
    for (std::size_t i = 0; i < spx.size(); ++i) {
        if (i != except_index)
            detail::tuple_add(result, spx[i]);
    }
    detail::tuple_divide_scalar(result, spx.size() - 1);
    return result;
}

template<class Fn, class Simplex>
typename Simplex::value_type nelder_mead(
        Fn f,
        const Simplex& initial_simplex,
        int max_iter,
        double ref_factor, double exp_factor, double con_factor)
{
    using std::begin;
    using std::end;

    auto trial_simplex(initial_simplex);
    std::array<typename std::result_of<Fn(typename Simplex::value_type)>::type,
        std::tuple_size<Simplex>::value> result;
    std::transform(begin(trial_simplex), end(trial_simplex), begin(result), f);

    auto extrema = std::minmax_element(begin(result), end(result));
    std::size_t best = std::distance(begin(result), extrema.first),
        worst = std::distance(begin(result), extrema.second);

    return trial_simplex[best];
}

}

#define CONVEX_HPP
#endif
