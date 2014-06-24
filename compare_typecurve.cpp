/*
 * Typecurve Aggregation and Fitting Technique Comparison
 *
 * dwt | terminus data science, LLC
 *
 * In this program, we will consider several techniques for aggregating
 * well production data reported on a monthly basis and fitting "type curves"
 * to represent the performance of the "average well" in the set.
 *
 * While we will work with Arps hyperbolic decline curves, the techniques
 * and conclusions here are generally applicable.
 */

#include <iostream>
#include <vector>
#include <cstddef>
#include <random>
#include <cmath>

#include "exponential.hpp"
#include "hyperbolic.hpp"
#include "hyptoexp.hpp"
#include "decline.hpp"
#include "bestfit.hpp"
#include "production.hpp"

const double year_days = 365.25;
const double month_days = 30.4;

const int forecast_years = 5;
const unsigned n_wells = 100;

int main()
{
    /*
     * We will begin by generating a random set of "true declines",
     * representing the (unknown) true performance of the wells in our group
     * of interest.
     *
     * We will give q_i a log-normal distribution with mean of 75 bbl/d and
     * log s.d. of 0.5 log bbl/d, D_i a normal distribution with mean 75% sec.
     * eff. / year and s.d. of 5% sec. eff. / year, and b a log-normal
     * distribution with mean 1.2 and log s.d. of 0.4.
     *
     * These parameters should generate the sorts of decline curves common
     * in tight liquids plays.
     */

    std::default_random_engine rand;
    std::lognormal_distribution<> qi_dist(std::log(75), 0.5);
    std::normal_distribution<> Di_dist(0.75, 0.05);
    std::lognormal_distribution<> b_dist(std::log(1.2), 0.4);

    std::vector<dca::arps_hyperbolic> true_declines;
    for (std::size_t i = 0; i < n_wells; ++i) {
        double b = b_dist(rand);
        true_declines.emplace_back(
                qi_dist(rand) /* bbl/d */ * year_days,
                dca::decline<dca::secant_effective>(Di_dist(rand), b),
                b);
    }

    for (const auto& decl: true_declines)
        std::cout << decl << '\n';

    /*
     * Now, let's generate three years worth of sample data using this decline.
     * We'll generate a series of times from t = 0 (beginning of first month) to
     * t = 2+11/12 (beginning of last month), as well as both a series of
     * *instantaneous rate* data, and a series of *interval volume* data.
     */

    std::vector<double> time(forecast_years * 12);
    dca::step_series(time.begin(), time.end(), 0.0, 1.0 / 12);
}
