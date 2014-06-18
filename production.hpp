#ifndef PRODUCTION_HPP

#include <cstddef>
#include <tuple>
#include <iterator>
#include <algorithm>
#include <utility>
#include <vector>

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

template<class AggFn, class ProdStreamIter, class OutIter>
inline OutIter aggregate_production(
        ProdStreamIter prod_begin, ProdStreamIter prod_end, OutIter out,
        std::size_t min_streams)
{
    using std::begin;
    using std::end;
    using prod_cont =
        typename std::iterator_traits<ProdStreamIter>::value_type;
    using prod_it = decltype(begin(std::declval<prod_cont>()));

    std::vector<std::pair<prod_it, prod_it>> streams;
    std::transform(prod_begin, prod_end, std::back_inserter(streams),
            [](const prod_cont& cont) {
                return std::make_pair(begin(cont), end(cont));
            });

    while (true) {
        std::size_t active_streams = 0;
        AggFn aggregate;

        for (auto& stream : streams)
            if (stream.first != stream.second) {
                ++active_streams;
                aggregate(*stream.first++);
            }

        if (active_streams < min_streams)
            break;

        *out = aggregate();
    }

    return out;
}

class mean {
    public:
        mean() noexcept : sum_(0.0), count_(0) {}

        void operator()(double x) noexcept
        {
            sum_ += x;
            ++count_;
        }

        double operator()() const noexcept
        {
            return sum_ / count_;
        }

    private:
        double sum_;
        std::size_t count_;
};

}

#define PRODUCTION_HPP
#endif
