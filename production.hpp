#ifndef PRODUCTION_HPP

#include <cstddef>
#include <tuple>
#include <iterator>
#include <algorithm>

namespace dca {

namespace detail {

template<std::size_t N>
struct shift_to_peak_impl {
    template<class TiedRateIter0, class Tuple, class... TiedRateIters>
    static void shift_to_peak(std::size_t shift, Tuple& result,
            TiedRateIter0 begin, TiedRateIters... rest)
    {
        std::get<std::tuple_size<Tuple>::value - N>(result) = begin + shift;
        shift_to_peak_impl<N - 1>::shift_to_peak(shift, result, rest...);
    }
};

template<>
struct shift_to_peak_impl<0> {
    template<class Tuple, class... TiedRateIters>
    static void shift_to_peak(std::size_t, Tuple&, TiedRateIters...) { }
};

}

template<class RateIter, class... TiedRateIters>
inline std::tuple<RateIter, TiedRateIters...>
shift_to_peak(RateIter major_begin, RateIter major_end,
        TiedRateIters... minor_begin)
{
    auto max = std::max_element(major_begin, major_end);
    auto shift = static_cast<std::size_t>(std::distance(major_begin, max));

    std::tuple<RateIter, TiedRateIters...> result;
    std::get<0>(result) = max;

    detail::shift_to_peak_impl<sizeof...(TiedRateIters)>::shift_to_peak(
            shift, result, minor_begin...);

    return result;
}

}

#define PRODUCTION_HPP
#endif
