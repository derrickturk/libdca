#include "dca/decline.hpp"

#define BOOST_TEST_MODULE decline
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

const double tolerance_pct = 1e-2;

BOOST_AUTO_TEST_CASE( conversions )
{
    using namespace dca;

    double dnom = 0.95, dtan = 0.15, dsec = 0.78, b = 1.5;

    // test round-trips from/to nominal
    BOOST_CHECK_CLOSE(dnom,
            (decline<secant_effective>(
                convert_decline<nominal, secant_effective>(dnom, b),
                b)),
            tolerance_pct);

    BOOST_CHECK_CLOSE(dnom,
            (decline<tangent_effective>(
                convert_decline<nominal, tangent_effective>(dnom)
                )),
            tolerance_pct);

    BOOST_CHECK_EQUAL(dnom, (decline<nominal>(dnom)));

    // test round-trips from/to tangent effective
    BOOST_CHECK_CLOSE(dtan,
            (convert_decline<secant_effective, tangent_effective>(
                convert_decline<tangent_effective, secant_effective>(
                    dtan, b), b)),
            tolerance_pct);

    BOOST_CHECK_CLOSE(dtan,
            (convert_decline<nominal, tangent_effective>(
                decline<tangent_effective>(dtan))),
            tolerance_pct);

    BOOST_CHECK_EQUAL(dtan,
            (convert_decline<tangent_effective, tangent_effective>(dtan)));

    // test round-trips from/to secant effective
    BOOST_CHECK_CLOSE(dsec,
            (convert_decline<tangent_effective, secant_effective>(
                convert_decline<secant_effective, tangent_effective>(
                    dsec, b), b)),
            tolerance_pct);

    BOOST_CHECK_CLOSE(dsec,
            (convert_decline<nominal, secant_effective>(
                decline<secant_effective>(dsec, b), b)),
            tolerance_pct);

    BOOST_CHECK_EQUAL(dsec,
            (convert_decline<secant_effective, secant_effective>(dsec, b)));
}
