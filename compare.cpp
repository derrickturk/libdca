/*
 * Decline Curve Best Fit and Typecurve Aggregation Technique Comparison
 */

#include <iostream>
#include <vector>
#include <cstddef>
#include <random>
#include <numeric>

#define DCA_IOSTREAMS

#include "exponential.hpp"
#include "hyperbolic.hpp"
#include "hyptoexp.hpp"
#include "decline.hpp"
#include "bestfit.hpp"
#include "production.hpp"

#undef DCA_IOSTREAMS

const double year_days = 365.25;
const double month_days = 30.4;

const int forecast_years = 5;

int main()
{
    /*
     * Let's begin by creating a "known" decline, which we'll use to generate
     * test data. We can then assess the ability of various techniques to
     * recover the true decline parameters given "fuzzed" data generated
     * using this decline.
     */

    auto true_decline = dca::arps_hyperbolic(
            /* qi = */ 100 * year_days, /* bbl/yr */
            /* Di = */ dca::decline<dca::secant_effective>(0.75 /* % / yr */,
                 1.5),
            /* b = */ 1.5);

    /*
     * Now, let's generate three years worth of sample data using this decline.
     * We'll generate a series of times from t = 0 (beginning of first month) to
     * t = 2+11/12 (beginning of last month), as well as both a series of
     * *instantaneous rate* data, and a series of *interval volume* data.
     */

    std::vector<double> time(forecast_years * 12);
    dca::step_series(time.begin(), time.end(), 0.0, 1.0 / 12);

    std::vector<double> true_instantaneous(time.size());
    std::transform(time.begin(), time.end(), true_instantaneous.begin(),
            [&](double t) { return true_decline.rate(t); });

    std::vector<double> true_interval(time.size());
    dca::interval_volumes(true_decline, true_interval.begin(), 0.0, 1.0 / 12,
            time.size());

    /* Let's take a look at the data. */
    std::cerr << "Data for true decline " << true_decline
        << "\nTime\tInstantaneous Rate (bbl/d)\tMonthly Volume (bbl)\n";
    for (std::size_t i = 0, sz = time.size(); i < sz; ++i)
        std::cerr << time[i] << '\t' << true_instantaneous[i] / year_days << '\t'
            << true_interval[i] << '\n';

    /*
     * In reality, things are more complicated. The peak rate will not likely
     * fall on the date of first production. Let's add a ramp-up period of
     * two months prior to the initiation of our hyperbolic decline.
     */
    true_instantaneous.insert(true_instantaneous.begin(), {
            25 /* bbl/d */ * year_days,
            50 /* bbl/d */ * year_days
        });
    true_instantaneous.resize(time.size());

    true_interval.insert(true_interval.begin(), {
            1000 /* bbl */,
            2200 /* bbl */
        });
    true_interval.resize(time.size());

    /* Now let's look at the data again. */
    std::cerr << "\nData for true decline " << true_decline << " with ramp-up"
        << "\nTime\tInstantaneous Rate (bbl/d)\tMonthly Volume (bbl)\n";
    for (std::size_t i = 0, sz = time.size(); i < sz; ++i)
        std::cerr << time[i] << '\t' << true_instantaneous[i] / year_days << '\t'
            << true_interval[i] << '\n';

    /*
     * When we fit this data, let's assume we only have access to the first
     * eight months (ramp-up and peak + 5) of producing history, and that
     * we only have access to the monthly data.
     *
     * Several approaches are possible. The most straightforward is to simply
     * find a decline fit which minimizes the error in the monthly interval
     * volumes. Because Arps declines are defined with the peak rate at time
     * t = 0, we need to begin our fit at the peak rate rather than at the
     * "actual" initial production date.
     */
    auto peak = std::get<0>(dca::shift_to_peak(
                true_interval.begin(), true_interval.end()));
    auto peak_shift = std::distance(true_interval.begin(), peak);
    auto fit_interval_with_shift =
        dca::best_from_interval_volume<dca::arps_hyperbolic>(
                peak, peak + 6 /* fit only peak and subsequent 5 months */,
                0.0, 1.0 / 12);

    /* Let's see what kind of fit we achieved. */
    std::cerr << "\nFit interval-volume peak+5 with shift-to-peak: "
        << fit_interval_with_shift << '\n';

    /*
     * Not bad! We recaptured the correct decline exactly.
     * Of course, real data is noisy, and the generative models don't exactly
     * match the models we are using to fit---we can take a look at this later.
     *
     * But for now, what about if we didn't shift to the peak before
     * trying to fit a model?
     */
    auto fit_interval_from_zero =
        dca::best_from_interval_volume<dca::arps_hyperbolic>(
                true_interval.begin(),
                /* fit ramp-up, peak, and subsequent 5 months */
                true_interval.begin() + 8,
                0.0, 1.0 / 12);

    /* How did we do? */
    std::cerr << "\nFit interval-volume out to peak+5 without shift-to-peak: "
        << fit_interval_from_zero << '\n';

    /*
     * That's not great news---we've underestimated the initial rate,
     * overstated the initial decline, overstated the b exponent, and generally
     * poorly fit the data.
     *
     * Another issue that often comes up has to do with the way available DCA
     * software treats monthly data.
     *
     * Or, more precisely, how it doesn't. Most software tools only operate on
     * instantaneous rate data, and assume any entered data are instantaneous
     * measurements.
     *
     * What happens if we let a tool fit a decline to our monthly
     * *interval volumes*, but assume that they reflect *instantaneous rates*?
     *
     * For this exercise, we'll pretend that the total production in a month
     * gets treated as an average which hits at the beginning of the month.
     * This reflects the way most DCA software works. We will grant ourselves
     * a shift to align the peak production at time 0.
     */
    std::vector<double> true_monthly_avg(true_interval.size());
    std::transform(true_interval.begin(), true_interval.end(),
            true_monthly_avg.begin(),
            [](double monthly) { return monthly / month_days * year_days; });
    auto fit_rate_from_average = dca::best_from_rate<dca::arps_hyperbolic>(
            true_monthly_avg.begin() + peak_shift,
            true_monthly_avg.begin() + peak_shift + 5,
            time.begin());

    /* How does that decline look? */
    std::cerr << "\nFit avg. monthly as instantaneous out to peak+5"
        " with shift-to-peak: " << fit_rate_from_average << '\n';

    /*
     * And what if we also neglect to shift everything to align to the peak rate?
     */
    auto fit_rate_from_average_no_shift = dca::best_from_rate<dca::arps_hyperbolic>(
            true_monthly_avg.begin(), true_monthly_avg.begin() + peak_shift + 5,
            time.begin());
    std::cerr << "\nFit avg. monthly as instantaneous out to peak+5"
        " without shift-to-peak: " << fit_rate_from_average_no_shift << '\n';

    /*
     * We can see that a big problem with fitting as if monthly volumes
     * were instantaneous rates is that we understate the true initial rate.
     * What if we adjust the fit's q_i to try to better match? In this case,
     * we'll assume that we know the true instantaneous peak rate from a
     * well test.
     */
    dca::arps_hyperbolic fit_rate_adjust_qi(
            true_instantaneous[peak_shift], /* use actual peak instantaneous */
            fit_rate_from_average.Di(), /* retain other parameters */
            fit_rate_from_average.b()
    );
    std::cerr << "\nFit avg. monthly as instantaneous out to peak+5"
        " with shift-to-peak and adjust qi: " << fit_rate_adjust_qi << '\n';

    /*
     * This looks a little better, but we are going to overshoot on
     * on ultimate recovery because we've taken a too-shallow curve and
     * lifted it to begin at a higher initial rate.
     *
     * One technique that might be incorrectly applied here is to
     * attempt to iterate on the b exponent in order to match some
     * "known" EUR value.
     *
     * In this case, we'll assume we know the true EUR, and adjust the
     * b exponent for our qi-adjusted curve to achieve the same EUR.
     */
    auto true_eur = dca::eur(true_decline,
            1.0 /* bbl/d */ * year_days,
            30 /* years */);
    auto new_b = std::get<0>(convex::nelder_mead([&](double b) {
            try {
                return std::pow(dca::eur(dca::arps_hyperbolic(
                            fit_rate_adjust_qi.qi(), fit_rate_adjust_qi.Di(), b),
                        1.0 * year_days, 30) - true_eur, 2);
            } catch (...) {
                return std::numeric_limits<double>::infinity();
            }
        }, convex::simplex<double> {
            { std::make_tuple(0.0), std::make_tuple(100) }
        }, 300));
    dca::arps_hyperbolic fit_rate_adjust_qi_b(
            fit_rate_adjust_qi.qi(),
            fit_rate_adjust_qi.Di(),
            new_b);
    std::cerr << "\nFit avg. monthly as instantaneous out to peak+5"
        " with shift-to-peak and adjust qi and b: "
        << fit_rate_adjust_qi_b << '\n';

    /*
     * Another "easy" way to compensate, as implemented by e.g. Fekete, is to
     * time-shift the monthly averages to the middle of the month, and fit
     * those values as instantaneous rates.
     */
    std::vector<double> fekete_time(time.size());
    std::transform(time.begin(), time.end(), fekete_time.begin(),
            [](double t) { return t + 1.0 / 24; });
    auto fit_fekete = dca::best_from_rate<dca::arps_hyperbolic>(
            true_monthly_avg.begin() + peak_shift,
            true_monthly_avg.begin() + peak_shift + 5,
            fekete_time.begin());
    std::cerr << "\nFit avg. monthly as instantaneous out to peak+5"
        " with shift-to-peak and \"Fekete time\": "
        << fit_fekete << '\n';

    /*
     * Let's compare these declines on an actual-vs-forecasted basis, out to
     * the end of the three year period.
     *
     * Regardless of how the curves were fit, we'll apply them beginning at
     * the peak rate, for a more fair comparison of how they'd actually
     * be used in forecasting.
     */

    std::cout << "\nTime\tCase\tType\tRate\tEUR\n";
    for (std::size_t i = 0, sz = time.size(); i < sz; ++i) {
        std::cout << time[i] << "\tActual\tIntervalAvg\t"
            << true_interval[i] / month_days << '\t' << true_eur << '\n';
        std::cout << time[i] << "\tActual\tInstantaneous\t"
            << true_instantaneous[i] / year_days << '\t' << true_eur << '\n';
    }

    std::vector<double> fit_interval(time.size());
    std::vector<std::pair<const char*, const dca::arps_hyperbolic*>> declines {
        { "IntervalFitShiftPeak", &fit_interval_with_shift },
        { "IntervalFitFromZero", &fit_interval_from_zero },
        { "RateFitShiftPeak", &fit_rate_from_average },
        { "RateFitFromZero", &fit_rate_from_average_no_shift },
        { "RateFitAdjustQi", &fit_rate_adjust_qi },
        { "RateFitAdjustQiB", &fit_rate_adjust_qi_b },
        { "RateFitFeketeTime", &fit_fekete }
    };

    for (const auto& decline : declines) {
        /* Leave ramp-up period at 0 */
        dca::interval_volumes(*decline.second,
                fit_interval.begin() + peak_shift, 0, 1.0 / 12,
                time.size() - peak_shift);

        auto eur = dca::eur(*decline.second, 1.0 * year_days, 30);

        for (std::size_t i = 0, sz = time.size(); i < sz; ++i) {
            std::cout << time[i] << '\t' << decline.first << "\tIntervalAvg\t"
                << fit_interval[i] / month_days << '\t' << eur << '\n';
            std::cout << time[i] << '\t' << decline.first << "\tInstantaneous\t"
                << decline.second->rate(time[i] - time[peak_shift]) / year_days
                << '\t' << eur << '\n';
        }
    }
}
