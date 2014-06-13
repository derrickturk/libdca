#include <iostream>
#include "exponential.hpp"
#include "hyperbolic.hpp"
#include "hyptoexp.hpp"

int main()
{
    dca::arps_exponential exp(1000,
            dca::decline<dca::tangent_effective>(0.95));

    std::cout << "exponential\n";
    for (double t = 0; t <= 12; t += 0.25)
        std::cout << "t = " << t << ", q = " << exp.rate(t) << ", Np = "
            << exp.cumulative(t) << '\n';

    dca::arps_hyperbolic hyp(1000,
            dca::decline<dca::tangent_effective>(0.95), 1.5);

    std::cout << "hyperbolic\n";
    for (double t = 0; t <= 12; t += 0.25)
        std::cout << "t = " << t << ", q = " << hyp.rate(t) << ", Np = "
            << hyp.cumulative(t) << ", D = " << hyp.D(t) << '\n';

    dca::arps_hyperbolic_to_exponential h2e(1000,
            dca::decline<dca::tangent_effective>(0.95), 1.5,
            dca::decline<dca::tangent_effective>(0.15));

    std::cout << "hyp2exp\n";
    for (double t = 0; t <= 12; t += 0.25)
        std::cout << "t = " << t << ", q = " << h2e.rate(t) << ", Np = "
            << h2e.cumulative(t) << ", D = " << h2e.D(t) << '\n';
}
