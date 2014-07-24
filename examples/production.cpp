#include "production.hpp"
#include <iterator>
#include <iostream>
#include <cassert>
#include <vector>

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

    std::vector<std::vector<double>> prod {
        { 4000, 3000, 2000, 1000, 500, 100 },
        { 1000,  750,  650,  500, 250,  50 },
        { 2500, 2000, 1250,  750, 500,  75 }
    };
    std::vector<double> type_well;

    dca::aggregate_production(prod.begin(), prod.end(),
            std::back_inserter(type_well), 3, dca::mean {});

    std::cout << "\nMean:\n";
    for (const auto& p : type_well)
        std::cout << p << '\n';

    dca::aggregate_production(prod.begin(), prod.end(),
            type_well.begin(), 3, dca::percentile(0.25));

    std::cout << "\nP25:\n";
    for (const auto& p : type_well)
        std::cout << p << '\n';

    std::vector<std::pair<double*, double*>> prod2 {
        std::make_pair(prod[0].data() + 1, prod[0].data() + 3),
        std::make_pair(prod[1].data() + 1, prod[1].data() + 3),
        std::make_pair(prod[2].data() + 1, prod[2].data() + 3)
    };

    type_well.clear();
    dca::aggregate_production(prod2.begin(), prod2.end(),
            std::back_inserter(type_well), 3, dca::mean {});

    std::cout << "\nMean:\n";
    for (const auto& p : type_well)
        std::cout << p << '\n';
}
