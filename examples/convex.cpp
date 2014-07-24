#include "convex.hpp"

#include <iostream>

inline double complex_fn(double x, double y)
{
    return -((x + 2) * (x + 2) - 3 * (y - 4) * (y - 4) * (y - 4) * (y - 4));
}

inline double t_complex_fn(const std::tuple<double, double>& t)
{
    return complex_fn(std::get<0>(t), std::get<1>(t));
}

int main()
{


    convex::simplex<double, double> s {
        {
            std::make_tuple(-2.5, -2.5),
            std::make_tuple(0, 0),
            std::make_tuple(2.5, 2.5)
        }
    };
    auto res = convex::nelder_mead(complex_fn, s, 1000);
    auto t_res = convex::nelder_mead(t_complex_fn, s, 1000);

    std::cout << "x = " << std::get<0>(res) << ", y = " << std::get<1>(res)
        << '\n';
    std::cout << "x = " << std::get<0>(t_res) << ", y = " << std::get<1>(t_res)
        << '\n';
}
