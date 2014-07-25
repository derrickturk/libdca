#ifndef TUPLE_TOOLS_HPP
#define TUPLE_TOOLS_HPP

#include <cstddef>
#include <type_traits>
#include <tuple>
#include <utility>

namespace tuple {

namespace detail {

template<std::size_t N>
struct apply_impl {
    template<class F, class Tuple, class... TParams>
    static auto apply(F&& fn, Tuple&& t, TParams&&... args)
      ->  decltype(
              apply_impl<N - 1>::apply(
                  std::forward<F>(fn), std::forward<Tuple>(t),
                  std::get<N - 1>(std::forward<Tuple>(t)),
                  std::forward<TParams>(args)...
          ))
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
    static auto apply(F&& fn, Tuple&&, TParams&&... args)
      -> decltype(std::forward<F>(fn)(std::forward<TParams>(args)...))
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

}

template<class F, class Tuple>
auto apply(F&& fn, Tuple&& t)
  -> decltype(detail::apply_impl<
          std::tuple_size<typename std::decay<Tuple>::type>::value
        >::apply(std::forward<F>(fn), std::forward<Tuple>(t)))
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

}

#endif
