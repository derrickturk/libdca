#ifndef CONVEX_HPP

#include <utility>
#include <type_traits>
#include <array>
#include <algorithm>
#include <cstddef>

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

}

template<class... Params>
using simplex =
    typename detail::simplex_traits<std::tuple<Params...>>::simplex_type;

template<class Simplex>
typename Simplex::value_type centroid(const Simplex& spx,
        std::size_t except_index);

template<class Fn, class Simplex>
typename Simplex::value_type nelder_mead(
        Fn f,
        const Simplex& initial_simplex,
        int max_iter,
        double ref_factor = 1.0,
        double exp_factor = 2.0,
        double con_factor = 0.5);

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
    auto cent = centroid(trial_simplex, worst);

    static_cast<void>(cent);
    static_cast<void>(max_iter);
    static_cast<void>(ref_factor);
    static_cast<void>(exp_factor);
    static_cast<void>(con_factor);

    return trial_simplex[best];
}

}

#define CONVEX_HPP
#endif
