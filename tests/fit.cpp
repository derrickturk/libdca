#include "dca/decline.hpp"
#include "dca/exponential.hpp"
#include "dca/hyperbolic.hpp"
#include "dca/hyptoexp.hpp"
#include "dca/bestfit.hpp"

#define BOOST_TEST_MODULE fit
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include <random>
#include <cmath>
#include <utility>
#include <vector>
#include <algorithm>

const double tolerance_pct = 1e-2;
const int n_test = 100;

template<class Decline>
auto forecast(const Decline& decline, double time_begin, double time_step,
        unsigned steps)
{
    std::vector<double> rate(steps), time(steps);
    double t = time_begin;
    std::generate(begin(time), end(time), [&]() {
        double current = t; t += time_step; return current;
    });
    std::transform(begin(time), end(time), begin(rate), [&](double t) {
        return decline.rate(t);
    });

    return make_pair(std::move(rate), std::move(time));
}

BOOST_AUTO_TEST_SUITE( fit_recovery )

BOOST_AUTO_TEST_CASE( exponential )
{
    std::mt19937 rng;

    std::uniform_real_distribution<> qi_log_dist(0.0, 7.0);
    std::uniform_real_distribution<> D_tangent_dist(0.0, 1.0);
    for (int i = 0; i < n_test; ++i) {
        dca::arps_exponential decl(std::pow(10.0, qi_log_dist(rng)),
                dca::decline<dca::tangent_effective>(D_tangent_dist(rng)));
        auto projection = forecast(decl, 0.0, 0.5, 100);
        auto fit = dca::best_from_rate<dca::arps_exponential>(
                begin(projection.first), end(projection.first),
                begin(projection.second));
        BOOST_CHECK_CLOSE(decl.qi(), fit.qi(), tolerance_pct);
        BOOST_CHECK_CLOSE(decl.D(), fit.D(), tolerance_pct);
    }
}

BOOST_AUTO_TEST_CASE( hyperbolic )
{
    std::mt19937 rng;

    std::uniform_real_distribution<> qi_log_dist(0.0, 7.0);
    std::uniform_real_distribution<> Di_tangent_dist(0.0, 1.0);
    std::uniform_real_distribution<> b_dist(0.0, 2.5);
    for (int i = 0; i < n_test; ++i) {
        dca::arps_hyperbolic decl(std::pow(10.0, qi_log_dist(rng)),
                dca::decline<dca::tangent_effective>(Di_tangent_dist(rng)),
                b_dist(rng));
        auto projection = forecast(decl, 0.0, 0.5, 100);
        auto fit = dca::best_from_rate<dca::arps_hyperbolic>(
                begin(projection.first), end(projection.first),
                begin(projection.second));
        BOOST_CHECK_CLOSE(decl.qi(), fit.qi(), tolerance_pct);
        BOOST_CHECK_CLOSE(decl.Di(), fit.Di(), tolerance_pct);
        BOOST_CHECK_CLOSE(decl.b(), fit.b(), tolerance_pct);
    }
}

BOOST_AUTO_TEST_CASE( hyptoexp )
{
    std::mt19937 rng;
}

BOOST_AUTO_TEST_SUITE_END()
