#ifndef EXPONENTIAL_HPP
#define EXPONENTIAL_HPP

#include <stdexcept>
#include <cmath>
#ifndef DCA_NO_IOSTREAMS
#include <iostream>
#endif

namespace dca {

class arps_exponential {
    public:
        arps_exponential(double qi, double D);

        const double& qi() const noexcept;
        const double& D() const noexcept;

        double rate(double time) const noexcept;
        double cumulative(double time) const noexcept;

    private:
        double qi_;
        double D_;

        static constexpr double eps_ = 1e-5;
};

inline arps_exponential::arps_exponential(double qi, double D)
    : qi_(qi), D_(D)
{
    if (qi_ < 0.0)
        throw std::out_of_range("qi must be non-negative.");
    if (D_ < 0.0)
        throw std::out_of_range("D must be non-negative.");
}

inline const double& arps_exponential::qi() const noexcept
{
    return qi_;
}

inline const double& arps_exponential::D() const noexcept
{
    return D_;
}

inline double arps_exponential::rate(double time) const noexcept
{
    if (time < 0.0) return 0.0;
    return qi_ * std::exp(-D_ * time);
}

inline double arps_exponential::cumulative(double time) const noexcept
{
    if (time < 0.0) return 0.0;
    if (D_ < eps_)
        return qi_ * time;
    return qi_ / D_ * (1.0 - std::exp(-D_ * time));
}

#ifndef DCA_NO_IOSTREAMS
inline std::ostream& operator<<(std::ostream& os, const arps_exponential& d)
{
    return os << "<Arps exponential decline: (qi = " << d.qi() << ", D = "
        << d.D() << ")>";
}
#endif

}

#endif
