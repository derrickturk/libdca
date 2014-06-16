#include "production.hpp"
#include <iterator>
#include <iostream>
#include <cassert>

int main()
{
    double oil[] = {
        1000,
        1200,
        950,
        750,
        500,
        100,
        50
    };

    double gas[] = {
        5000,
        4000,
        3000,
        2000,
        1000,
        500,
        250
    };

    using std::begin;
    using std::end;
    auto from_peak = dca::shift_to_peak(begin(oil), end(oil),
            begin(gas));

    for (auto o = std::get<0>(from_peak), g = std::get<1>(from_peak);
            o != end(oil); ++o, ++g)
        std::cout << "Oil: " << *o << ", Gas: " << *g << '\n';

    auto oil_from_peak = dca::shift_to_peak(begin(oil), end(oil));
    assert(std::get<0>(oil_from_peak) == std::get<0>(from_peak));

    auto multi = dca::shift_to_peak(begin(oil), end(oil),
            begin(gas), begin(gas), begin(gas));
    assert(std::get<3>(multi) == std::get<1>(multi));
}
