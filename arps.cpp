#include <iostream>
#include "exponential.hpp"
#include "hyperbolic.hpp"

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
            << hyp.cumulative(t) << '\n';
}
