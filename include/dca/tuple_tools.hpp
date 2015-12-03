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

}

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

template<class>
struct result_of;

template<class F, class Tuple>
struct result_of<F(Tuple)>
{
    typedef decltype(apply(std::declval<F>(), std::declval<Tuple>())) type;
};

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

}

#endif
