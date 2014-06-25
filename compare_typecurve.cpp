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
#include <iterator>
#include <numeric>

#include "exponential.hpp"
#include "hyperbolic.hpp"
#include "hyptoexp.hpp"
#include "decline.hpp"
#include "bestfit.hpp"
#include "production.hpp"

const double year_days = 365.25;
const double month_days = 30.4;

const int forecast_years = 5;
const int fit_months = 6;
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

    /*
     * Now, let's generate a suitable amount of (post-peak) data for each well.
     */

    std::vector<double> time(forecast_years * 12);
    dca::step_series(time.begin(), time.end(), 0.0, 1.0 / 12);

    std::vector<std::vector<double>> rate(true_declines.size());
    std::transform(true_declines.begin(), true_declines.end(), rate.begin(),
            [&](const dca::arps_hyperbolic& decl) -> std::vector<double> {
                std::vector<double> result(forecast_years * 12);
                std::transform(time.begin(), time.end(), result.begin(),
                    [&](double t) { return decl.rate(t); });
                return result;
            });

    std::vector<std::vector<double>> production(true_declines.size());
    std::transform(true_declines.begin(), true_declines.end(),
            production.begin(),
            [&](const dca::arps_hyperbolic& decl) -> std::vector<double> {
                std::vector<double> result(forecast_years * 12);
                dca::interval_volumes(decl, result.begin(), 0.0, 1.0 / 12,
                    forecast_years * 12);
                return result;
            });

    /*
     * We'll also calculate EURs and 3-year cumulatives per well, for future
     * comparison to typecurve results.
     */
    std::vector<double> true_eur(true_declines.size()),
        true_3_cum(true_declines.size());
    std::transform(true_declines.begin(), true_declines.end(),
            true_eur.begin(), [](const dca::arps_hyperbolic& decl) {
                return dca::eur(decl, 1.0 * year_days, 30);
            });
    std::transform(true_declines.begin(), true_declines.end(),
            true_3_cum.begin(), [](const dca::arps_hyperbolic& decl) {
                return decl.cumulative(3.0);
            });
    double true_avg_eur = std::accumulate(true_eur.begin(), true_eur.end(),
            0.0) / true_eur.size(),
           true_avg_3_cum = std::accumulate(true_3_cum.begin(),
                   true_3_cum.end(), 0.0) / true_3_cum.size();

    std::cerr << "Avg of EUR = " << true_avg_eur / 1000 << " Mbbl.\n"
        << "Avg of 3-year cum. = " << true_avg_3_cum / 1000 << " Mbbl.\n";

    /*
     * Having previously established the correctness of the technique, we'll
     * apply an interval-volume shift-to-peak hyperbolic fit to each well,
     * using only the first six months of data.
     */

    std::vector<dca::arps_hyperbolic> fit_declines;
    std::vector<double> interval(6);
    std::transform(production.begin(), production.end(),
            std::back_inserter(fit_declines),
            [&](const std::vector<double>& prod) {
                return dca::best_from_interval_volume<dca::arps_hyperbolic>(
                    prod.begin(), prod.begin() + fit_months, 0.0, 1.0 / 12);
            });

    /*
     * We might attempt to use the average of each decline parameter in order
     * to generate an "average decline."
     */

    double avg_qi = 0.0, avg_Di = 0.0, avg_b = 0.0;
    for (const auto& decl : fit_declines) {
        avg_qi += decl.qi();
        avg_Di += decl.Di();
        avg_b += decl.b();
    }
    avg_qi /= fit_declines.size();
    avg_Di /= fit_declines.size();
    avg_b /= fit_declines.size();

    dca::arps_hyperbolic avg_params_decline(avg_qi, avg_Di, avg_b);

    std::cerr << "Avg params decline: " << avg_params_decline << '\n';
    std::cerr << "EUR = "
        << dca::eur(avg_params_decline, 1.0 * year_days, 30) / 1000
        << " Mbbl.\n" << "3-year cum. = "
        << avg_params_decline.cumulative(3) / 1000 << " Mbbl.\n";

    /*
     * That doesn't seem to match our known distributions very well. The better
     * approach is to aggregate the average production across wells, then fit
     * the desired decline curve model to the aggregate production.
     */

    std::vector<double> avg_rate;
    dca::aggregate_production(rate.begin(), rate.end(),
            std::back_inserter(avg_rate),
            rate.size() / 2, /* until less than half of wells producing */
            dca::mean {});

    std::vector<double> avg_production;
    dca::aggregate_production(production.begin(), production.end(),
            std::back_inserter(avg_production),
            production.size() / 2, /* until less than half of wells producing */
            dca::mean {});

    auto avg_prod_decline =
        dca::best_from_interval_volume<dca::arps_hyperbolic>(
            avg_production.begin(), avg_production.begin() + fit_months,
            0.0, 1.0 / 12);

    std::cerr << "Avg production decline: " << avg_prod_decline << '\n';
    std::cerr << "EUR = "
        << dca::eur(avg_prod_decline, 1.0 * year_days, 30) / 1000
        << " Mbbl.\n" << "3-year cum. = "
        << avg_prod_decline.cumulative(3) / 1000 << " Mbbl.\n";

    /*
     * As you can see, there is a significant difference between the curves.
     * Let's take a look.
     */

    std::cout << "Time\tCase\tType\tRate\n";
    for (std::size_t i = 0, sz = time.size(); i < sz; ++i) {
        std::cout << time[i] << "\tActualAvg\tIntervalAvg\t"
            << avg_production[i] / month_days << '\n';
        std::cout << time[i] << "\tActualAvg\tInstantaneous\t"
            << avg_rate[i] / year_days << '\n';
    }

    std::vector<double> fit_interval(time.size());
    std::vector<std::pair<const char*, const dca::arps_hyperbolic*>> declines {
        { "AvgParams", &avg_params_decline },
        { "AvgProd", &avg_prod_decline }
    };

    for (const auto& decline : declines) {
        /* Leave ramp-up period at 0 */
        dca::interval_volumes(*decline.second,
                fit_interval.begin(), 0, 1.0 / 12, time.size());

        for (std::size_t i = 0, sz = time.size(); i < sz; ++i) {
            std::cout << time[i] << '\t' << decline.first << "\tIntervalAvg\t"
                << fit_interval[i] / month_days << '\n';
            std::cout << time[i] << '\t' << decline.first << "\tInstantaneous\t"
                << decline.second->rate(time[i]) / year_days << '\n';
        }
    }
}
