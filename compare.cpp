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

int main()
{
    /*
     * Let's begin by creating a "known" decline, which we'll use to generate
     * test data. We can then assess the ability of various techniques to
     * recover the true decline parameters given "fuzzed" data generated
     * using this decline.
     */

    auto true_decline = dca::arps_hyperbolic(
            /* qi = */ 100 * 365.25, /* bbl/yr */
            /* Di = */ dca::decline<dca::secant_effective>(0.75 /* % / yr */,
                 1.5),
            /* b = */ 1.5);

    /*
     * Now, let's generate three years worth of sample data using this decline.
     * We'll generate a series of times from t = 0 (beginning of first month) to
     * t = 2+11/12 (beginning of last month), as well as both a series of
     * *instantaneous rate* data, and a series of *interval volume* data.
     */

    std::vector<double> time(3 * 12);
    dca::step_series(time.begin(), time.end(), 0.0, 1.0 / 12);

    std::vector<double> true_instantaneous(time.size());
    std::transform(time.begin(), time.end(), true_instantaneous.begin(),
            [&](double t) { return true_decline.rate(t); });

    std::vector<double> true_interval(time.size());
    dca::interval_volumes(true_decline, true_interval.begin(), 0.0, 1.0 / 12,
            time.size());

    /* Let's take a look at the data. */
    std::cout << "Data for true decline " << true_decline
        << "\nTime\tInstantaneous Rate (bbl/d)\tMonthly Volume (bbl)\n";
    for (std::size_t i = 0, sz = time.size(); i < sz; ++i)
        std::cout << time[i] << '\t' << true_instantaneous[i] / 365.25 << '\t'
            << true_interval[i] << '\n';

    /*
     * In reality, things are more complicated. The peak rate will not likely
     * fall on the date of first production. Let's add a ramp-up period of
     * two months prior to the initiation of our hyperbolic decline.
     */

    time.resize(8);
    dca::step_series(time.begin(), time.end(), 0.0, 1.0 / 12);

    true_instantaneous.insert(true_instantaneous.begin(), {
            25 /* bbl/d */ * 365.25,
            50 /* bbl/d */ * 365.25
        });

    true_interval.insert(true_interval.begin(), {
            1000 /* bbl */,
            2200 /* bbl */
        });

    /* Now let's look at the data again. */
    std::cout << "\nData for true decline " << true_decline << " with ramp-up"
        << "\nTime\tInstantaneous Rate (bbl/d)\tMonthly Volume (bbl)\n";
    for (std::size_t i = 0, sz = time.size(); i < sz; ++i)
        std::cout << time[i] << '\t' << true_instantaneous[i] / 365.25 << '\t'
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
    std::cout << "\nFit interval-volume peak+5 with shift-to-peak: "
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
    std::cout << "\nFit interval-volume out to peak+5 without shift-to-peak: "
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
            [](double monthly) { return monthly / 30.4 * 365.25; });
    auto fit_rate_from_average = dca::best_from_rate<dca::arps_hyperbolic>(
            true_monthly_avg.begin() + peak_shift,
            true_monthly_avg.begin() + peak_shift + 5,
            time.begin() + peak_shift);

    /* How does that decline look? */
    std::cout << "\nFit avg. monthly as instantaneous out to peak+5"
        " with shift-to-peak: " << fit_rate_from_average << '\n';

    /*
     * Let's compare these declines on an actual-vs-forecasted basis, out to
     * the end of the three year period.
     */

    std::cout << "Time\tCase\tType\tRate\n";
    for (std::size_t i = 0, sz = time.size(); i < sz; ++i) {
        std::cout << time[i] << "\tActual\tIntervalAvg\t"
            << true_interval[i] / 30.4 << '\n';
        std::cout << time[i] << "\tActual\tInstantaneous\t"
            << true_instantaneous[i] / 365.25 << '\n';
    }

    std::vector<double> fit_interval(time.size());

    /* Leave ramp-up period at 0 */
    dca::interval_volumes(fit_interval_with_shift,
            fit_interval.begin() + peak_shift,
            0, 1.0 / 12, time.size() - peak_shift);
    for (std::size_t i = 0, sz = time.size(); i < sz; ++i) {
        std::cout << time[i] << "\tIntervalFitShiftPeak\tIntervalAvg\t"
            << fit_interval[i] / 30.4 << '\n';
        std::cout << time[i] << "\tIntervalFitShiftPeak\tInstantaneous\t"
            << fit_interval_with_shift.rate(time[i] - time[peak_shift]) / 365.25
            << '\n';
    }

    dca::interval_volumes(fit_interval_from_zero,
            fit_interval.begin(), 0, 1.0 / 12, time.size());
    for (std::size_t i = 0, sz = time.size(); i < sz; ++i) {
        std::cout << time[i] << "\tIntervalFitFromZero\tIntervalAvg\t"
            << fit_interval[i] / 30.4 << '\n';
        std::cout << time[i] << "\tIntervalFitFromZero\tInstantaneous\t"
            << fit_interval_from_zero.rate(time[i]) / 365.25
            << '\n';
    }

    dca::interval_volumes(fit_rate_from_average,
            fit_interval.begin() + peak_shift, 0, 1.0 / 12,
            time.size() - peak_shift);
    for (std::size_t i = 0, sz = time.size(); i < sz; ++i) {
        std::cout << time[i] << "\tRateFitShiftPeak\tIntervalAvg\t"
            << fit_interval[i] / 30.4 << '\n';
        std::cout << time[i] << "\tRateFitShiftPeak\tInstantaneous\t"
            << fit_rate_from_average.rate(time[i] - time[peak_shift]) / 365.25
            << '\n';
    }

}
