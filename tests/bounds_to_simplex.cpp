#include "dca/convex.hpp"

#define BOOST_TEST_MODULE bounds_to_simplex
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/test_tools.hpp>

#include <utility>
#include <tuple>

#include <iostream>

/*
// WTF BOOST
namespace boost {
namespace test_tools {

template<class... Params>
struct print_log_value<std::tuple<Params...>> {
    void operator()(std::ostream& os, const std::tuple<Params...>& tuple)
    {
        using tuple::operator<<;
        os << tuple;
    }
};

}
}
*/

BOOST_AUTO_TEST_SUITE( inner_simplex )

BOOST_AUTO_TEST_CASE( dim_2 )
{
    using tuple::operator<<;
    auto is = convex::inner_simplex(std::make_pair(
                std::make_tuple(0.0, 0.0),
                std::make_tuple(10.0, 10.0)));
    BOOST_CHECK(is[0] == std::make_tuple(0.0, 0.0));
    BOOST_CHECK(is[1] == std::make_tuple(10.0, 0.0));
    BOOST_CHECK(is[2] == std::make_tuple(5.0, 10.0));
}

BOOST_AUTO_TEST_CASE( dim_3 )
{
    auto is = convex::inner_simplex(std::make_pair(
                std::make_tuple(0.0, 0.0, 0.0),
                std::make_tuple(10.0, 10.0, 10.0)));
    BOOST_CHECK(is[0] == std::make_tuple(0.0, 0.0, 0.0));
    BOOST_CHECK(is[1] == std::make_tuple(10.0, 0.0, 0.0));
    BOOST_CHECK(is[2] == std::make_tuple(5.0, 10.0, 0.0));
    BOOST_CHECK(is[3] == std::make_tuple(5.0, 5.0, 10.0));
}

BOOST_AUTO_TEST_CASE( dim_4 )
{
    auto is = convex::inner_simplex(std::make_pair(
                std::make_tuple(0.0, 0.0, 0.0, 0.0),
                std::make_tuple(10.0, 10.0, 10.0, 10.0)));
    BOOST_CHECK(is[0] == std::make_tuple(0.0, 0.0, 0.0, 0.0));
    BOOST_CHECK(is[1] == std::make_tuple(10.0, 0.0, 0.0, 0.0));
    BOOST_CHECK(is[2] == std::make_tuple(5.0, 10.0, 0.0, 0.0));
    BOOST_CHECK(is[3] == std::make_tuple(5.0, 5.0, 10.0, 0.0));
    BOOST_CHECK(is[4] == std::make_tuple(5.0, 5.0, 5.0, 10.0));
}

BOOST_AUTO_TEST_SUITE_END()
