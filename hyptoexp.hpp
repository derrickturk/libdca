#ifndef HYP2EXP_HPP

#include "decline.hpp"
#include "exponential.hpp"
#include "hyperbolic.hpp"
#include <stdexcept>
#include <cmath>

namespace dca {

class arps_hyperbolic_to_exponential :
  private arps_hyperbolic, private arps_exponential {
    public:
        arps_hyperbolic_to_exponential
            (double qi, double Di, double b, double Df);

        const double& qi() const noexcept;
        const double& Di() const noexcept;
        const double& b() const noexcept;
        const double& Df() const noexcept;

        double rate(double time) const noexcept;
        double cumulative(double time) const noexcept;
        double D(double time) const noexcept;

    private:
        double t_trans_;
};

inline arps_hyperbolic_to_exponential::arps_hyperbolic_to_exponential(
        double qi, double Di, double b, double Df)
    : arps_hyperbolic(qi, Di, b),
      arps_exponential(arps_hyperbolic::rate((Di / Df - 1.0) / (b * Di)), Df),
      t_trans_((Di / Df - 1.0) / (b * Di))
{
    if (Df <= 0) throw std::out_of_range("Df must be non-negative.");
    // note: if Df > Di, transition will occur at t < 0 and the curve
    // will be treated as wholly exponential
}

}

#define HYP2EXP_HPP
#endif
