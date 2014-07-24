#include "dca/tuple_tools.hpp"

#define BOOST_TEST_MODULE tuple_tools
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_CASE( tuple_apply )
{
    auto f = [](int x, int y) { return x + y; };
    auto g = [](int x) { return x; };
    auto h = []() { return 1; };

    BOOST_CHECK(tuple::apply(f, std::make_tuple(1, 2)) == 3);
    BOOST_CHECK(tuple::apply(g, std::make_tuple(1)) == 1);
    BOOST_CHECK(tuple::apply(h, std::make_tuple()) == 1);
}

BOOST_AUTO_TEST_CASE( tuple_construct )
{
    struct X {
        int x, y;
        X(int a, b) : x(a), b(b) { }
    };

    struct Y {
        int x;
        Y(int a) : x(a) { }
    };

    struct Z {
        Z() { }
    };

    auto x = tuple::construct<X>(std::make_tuple(1, 2));
    auto y = tuple::construct<Y>(std::make_tuple(1));
    auto z = tuple::construct<Z>(std::make_tuple());

    BOOST_CHECK(x.x == 1 && x.y == 2);
    BOOST_CHECK(y.x == 1);
}
