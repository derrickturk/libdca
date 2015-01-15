#include <iostream>
#include <vector>
#include "dca/decline.hpp"
#include "dca/exponential.hpp"
#include "dca/hyperbolic.hpp"
#include "dca/hyptoexp.hpp"
#include "dca/bestfit.hpp"
#include "dca/any_decline.hpp"

int main()
{
    dca::arps_exponential exp(1000,
            dca::decline<dca::tangent_effective>(0.95));

    std::cout << "exponential\n";
    for (double t = 0; t <= 12; t += 0.5)
        std::cout << "t = " << t << ", q = " << exp.rate(t) << ", Np = "
            << exp.cumulative(t) << '\n';
    std::cout << "EUR: " << dca::eur(exp, 1, 30) << '\n';

    dca::arps_hyperbolic hyp(1000,
            dca::decline<dca::tangent_effective>(0.95), 1.5);

    std::cout << "hyperbolic\n";
    for (double t = 0; t <= 12; t += 0.5)
        std::cout << "t = " << t << ", q = " << hyp.rate(t) << ", Np = "
            << hyp.cumulative(t) << ", D = " << hyp.D(t) << '\n';
    std::cout << "EUR: " << dca::eur(hyp, 1, 30) << '\n';

    dca::arps_hyperbolic_to_exponential h2e(1000,
            dca::decline<dca::tangent_effective>(0.95), 1.5,
            dca::decline<dca::tangent_effective>(0.15));

    std::cout << "hyp2exp\n";
    for (double t = 0; t <= 12; t += 0.5)
        std::cout << "t = " << t << ", q = " << h2e.rate(t) << ", Np = "
            << h2e.cumulative(t) << ", D = " << h2e.D(t) << '\n';
    std::cout << "EUR: " << dca::eur(h2e, 1, 30) << '\n';

    std::vector<double> time;
    for (double t = 0; t <= 25; t += 0.5)
        time.push_back(t);

    std::vector<double> rate, interval;
    for (auto t: time) {
        rate.push_back(h2e.rate(t));
        interval.push_back(h2e.cumulative(t + 0.5));
    }

    for (auto i = time.size() - 1; i > 0; --i)
        interval[i] -= interval[i - 1];

    auto best_rate = dca::best_from_rate<dca::arps_hyperbolic_to_exponential>(
            rate.begin(), rate.end(), time.begin());

    auto best_interval =
        dca::best_from_interval_volume<dca::arps_hyperbolic_to_exponential>(
                interval.begin(), interval.end(), 0, 0.5);

    std::cout << "best rate fit (" << best_rate.qi() << ", " << best_rate.Di()
        << ", " << best_rate.b() << ", " << best_rate.Df() << ")\n";
    std::cout << "best interval fit (" << best_interval.qi() << ", "
        << best_interval.Di() << ", " << best_interval.b() << ", "
        << best_interval.Df() << ")\n";

    dca::any decline(exp), decline2(hyp);
    decline = h2e;
    decline = decline2;

    std::cout << "any: " << decline << '\n';
    for (double t = 0; t <= 12; t += 0.5)
        std::cout << "t = " << t << ", q = " << decline.rate(t) << ", Np = "
            << decline.cumulative(t) << '\n';
    std::cout << "EUR: " << dca::eur(decline, 1, 30) << '\n';

    std::vector<double> intervals;
    dca::interval_volumes(decline, std::back_inserter(intervals),
            0, 1.0 / 12, 24);
    std::cout << "Interval volumes:\n";
    for (auto d: intervals)
        std::cout << d << '\n';
}
