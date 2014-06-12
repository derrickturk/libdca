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
    static const std::size_t simplex_length =
        1 + std::tuple_size<Params...>::value;
    using simplex_type = std::array<std::tuple<Params...>, simplex_length>;
};

}

template<class... Params>
using simplex = typename detail::simplex_traits<Params...>::simplex_type;

template<class... Params>
std::tuple<Params...> centroid(const simplex<Params...>& spx,
        std::size_t except_index);

template<class Fn, class... Params>
std::tuple<Params...> nelder_mead(
        Fn f,
        int max_iter,
        const simplex<Params...>& initial_simplex,
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

template<class Fn, class... Params>
std::tuple<Params...> nelder_mead(
        Fn f,
        const simplex<Params...>& initial_simplex,
        int max_iter,
        double ref_factor, double exp_factor, double con_factor)
{
    auto trial_simplex { initial_simplex };
    std::array<typename std::result_of<Fn(std::tuple<Params...>)>::type,
        simplex_traits<Params...>::simplex_length> result;
    std::transform(begin(simplex), end(simplex), begin(result), f);

    auto extrema = std::minmax_element(begin(result), end(result));
    std::size_t best = std::distance(begin(result), extrema.first),
        worst = std::distance(begin(result), extrema.second);
}

}

#define CONVEX_HPP
#endif
