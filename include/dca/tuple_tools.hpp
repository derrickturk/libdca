#ifndef TUPLE_TOOLS_HPP
#define TUPLE_TOOLS_HPP

#include <cstddef>
#include <type_traits>
#include <tuple>
#include <utility>
#ifndef DCA_NO_IOSTREAMS
#include <iostream>
#endif

namespace tuple {

namespace detail {

template<std::size_t N>
struct apply_impl {
    template<class F, class Tuple, class... TParams>
    static decltype(auto) apply(F&& fn, Tuple&& t, TParams&&... args)
    {
        return apply_impl<N - 1>::apply(
                std::forward<F>(fn), std::forward<Tuple>(t),
                std::get<N - 1>(std::forward<Tuple>(t)),
                std::forward<TParams>(args)...
                );
    }
};

template<>
struct apply_impl<0> {
    template<class F, class Tuple, class... TParams>
    static decltype(auto) apply(F&& fn, Tuple&&, TParams&&... args)
    {
        return std::forward<F>(fn)(std::forward<TParams>(args)...);
    }
};

template<std::size_t N>
struct construct_impl {
    template<class T, class Tuple, class... TParams>
    static T construct(Tuple&& t, TParams&&... args)
    {
        return construct_impl<N - 1>::template construct<T>(
                std::forward<Tuple>(t),
                std::get<N - 1>(std::forward<Tuple>(t)),
                std::forward<TParams>(args)...
                );
    }
};

template<>
struct construct_impl<0> {
    template<class T, class Tuple, class... TParams>
    static T construct(Tuple&&, TParams&&... args)
    {
        return T(std::forward<TParams>(args)...);
    }
};

template<class Tuple, std::size_t... I>
auto shuffle_left_construct(const Tuple& t, std::index_sequence<I...>)
{
    return std::make_tuple(std::get<I + 1>(t)..., std::get<0>(t));
}

template<class Tuple, std::size_t... I>
auto shuffle_right_construct(const Tuple& t, std::index_sequence<I...>)
{
    return std::make_tuple(std::get<std::tuple_size<Tuple>::value - 1>(t),
            std::get<I>(t)...);
}

template<class First, class... Params>
auto shuffle_left_impl(const std::tuple<First, Params...>& tuple)
{
    return shuffle_left_construct(tuple,
            std::make_index_sequence<sizeof...(Params)>());
}

template<class First, class... Params>
auto shuffle_right_impl(const std::tuple<First, Params...>& tuple)
{
    return shuffle_right_construct(tuple,
            std::make_index_sequence<sizeof...(Params)>());
}

#ifndef DCA_NO_IOSTREAMS
template<class Tuple, char Delim, std::size_t I>
struct streamto_impl {
    static std::ostream& streamto(std::ostream& os, const Tuple& tuple)
    {
        streamto_impl<Tuple, Delim, I - 1>::streamto(os, tuple);
        return os << Delim << std::get<I - 1>(tuple);
    }
};

template<class Tuple, char Delim>
struct streamto_impl<Tuple, Delim, 1>
{
    static std::ostream& streamto(std::ostream& os, const Tuple& tuple)
    {
        return os << std::get<0>(tuple);
    }
};

template<class Tuple, char Delim>
struct streamto_impl<Tuple, Delim, 0>
{
    static std::ostream& streamto(std::ostream& os, const Tuple& tuple)
    {
        return os;
    }
};
#endif

} // namespace detail

template<class F, class Tuple>
decltype(auto) apply(F&& fn, Tuple&& t)
{
    return detail::apply_impl<
        std::tuple_size<typename std::decay<Tuple>::type>::value
      >::apply(std::forward<F>(fn), std::forward<Tuple>(t));
}

template<class T, class Tuple>
T construct(Tuple&& t)
{
    return detail::construct_impl<
        std::tuple_size<typename std::decay<Tuple>::type>::value
      >::template construct<T>(std::forward<Tuple>(t));
}

template<class... Params>
auto shuffle_left(const std::tuple<Params...>& tuple)
{
    return detail::shuffle_left_impl(tuple);
}

template<>
auto shuffle_left(const std::tuple<>& tuple)
{
    return tuple;
}

template<class... Params>
auto shuffle_right(const std::tuple<Params...>& tuple)
{
    return detail::shuffle_right_impl(tuple);
}

template<>
auto shuffle_right(const std::tuple<>& tuple)
{
    return tuple;
}

#ifndef DCA_NO_IOSTREAMS
template<class... Params, char Begin='(', char End=')', char Delim=','>
std::ostream& operator<<(std::ostream& os, const std::tuple<Params...>& tuple)
{
    os << Begin;
    detail::streamto_impl<decltype(tuple), Delim, sizeof...(Params)>::
        streamto(os, tuple);
    return os << End;
}
#endif

template<class>
struct result_of;

template<class F, class Tuple>
struct result_of<F(Tuple)>
{
    typedef decltype(apply(std::declval<F>(), std::declval<Tuple>())) type;
};

template<class F, class Tuple>
using result_of_t = typename result_of<F(Tuple)>::type;

}

#endif
