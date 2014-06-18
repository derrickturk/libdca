#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <iterator>
#include <unordered_map>
#include <cstddef>
#include <algorithm>
#include <cstdlib>

#include "exponential.hpp"
#include "hyperbolic.hpp"
#include "hyptoexp.hpp"
#include "bestfit.hpp"
#include "production.hpp"

namespace params {

const std::string id_field = "PROPNUM";
static const std::string oil_field = "OIL";
static const std::string gas_field = "GAS";
static const auto aggregation = dca::mean {};

}

using dataset = std::unordered_map<std::string, std::vector<std::string>>;

dataset read_delimited(std::istream& is, char delim = '\t');

template<class F>
F foreach_well(const dataset& data, F fn,
        std::string id_field = params::id_field);

template<class I, class T, class F>
void for_delimited(I begin, I end, T delim, F fn)
{
    I pos;
    while ((pos = std::find(begin, end, delim)) != end) {
        fn(begin, pos);
        pos += 1;
        begin = pos;
    }
    if (begin != end)
        fn(begin, end);
}

int main(int argc, char* argv[])
{
    std::vector<std::string> args(argv, argv + argc);
    dataset data;

    if (args.size() == 2) {
        if (std::ifstream in { args[1] })
            data = read_delimited(in);
        else {
            std::cerr << "Unable to read from " << args[1] << '\n';
            return 0;
        }
    } else if (args.size() > 2) {
        std::cerr << "Usage: " << (args.empty() ? "bbep_fit" : args[0])
            << "<delim-file>\n";
    } else {
        data = read_delimited(std::cin);
    }

    std::vector<std::vector<double>> oil;
    std::vector<std::vector<double>> gas;
    foreach_well(data, [&](const dataset& well) {
        auto&& oil_text = well.at(params::oil_field);
        oil.emplace_back(oil_text.size());
        std::transform(oil_text.begin(), oil_text.end(), oil.back().begin(),
            [](const std::string& s) { return std::strtod(s.c_str(), nullptr); }
        );

        auto&& gas_text = well.at(params::gas_field);
        gas.emplace_back(gas_text.size());
        std::transform(gas_text.begin(), gas_text.end(), gas.back().begin(),
            [](const std::string& s) { return std::strtod(s.c_str(), nullptr); }
        );
    });

    double oil_avg_shift = 0;
    for (auto& oil_rec : oil) {
        auto peak = dca::shift_to_peak(oil_rec.begin(), oil_rec.end());
        oil_avg_shift += std::distance(oil_rec.begin(), std::get<0>(peak));
        oil_rec.erase(oil_rec.begin(), std::get<0>(peak));
    }
    oil_avg_shift /= oil.size();

    double gas_avg_shift = 0;
    for (auto& gas_rec : gas) {
        auto peak = dca::shift_to_peak(gas_rec.begin(), gas_rec.end());
        gas_avg_shift += std::distance(gas_rec.begin(), std::get<0>(peak));
        gas_rec.erase(gas_rec.begin(), std::get<0>(peak));
    }
    gas_avg_shift /= gas.size();

    std::vector<double> oil_tw, gas_tw;
    dca::aggregate_production(oil.begin(), oil.end(),
            std::back_inserter(oil_tw), std::floor(oil.size() / 3),
            params::aggregation);
    dca::aggregate_production(gas.begin(), gas.end(),
            std::back_inserter(gas_tw), std::floor(gas.size() / 3),
            params::aggregation);

    auto oil_tc = dca::best_from_interval_volume<dca::arps_hyperbolic>(
            oil_tw.begin(), oil_tw.end(), 0, 1.0 / 12.0);
    auto gas_tc = dca::best_from_interval_volume<dca::arps_hyperbolic>(
            gas_tw.begin(), gas_tw.end(), 0, 1.0 / 12.0);

    double t_eur;
    double oil_eur = dca::eur(
            dca::arps_hyperbolic_to_exponential(
                oil_tc.qi(),
                oil_tc.Di(),
                oil_tc.b(),
                dca::decline<dca::tangent_effective>(0.05)
            ),
            365.25, // bbl/yr = 1 bbl/day
            30, // years,
            &t_eur
    );

    double gas_eur = dca::arps_hyperbolic_to_exponential(
            gas_tc.qi(),
            gas_tc.Di(),
            gas_tc.b(),
            dca::decline<dca::tangent_effective>(0.05)
            ).cumulative(t_eur - (gas_avg_shift - oil_avg_shift));

    std::cout << "Oil Avg. Shift: " << oil_avg_shift << " months\n";
    std::cout << "Oil Type Well:\nMonth\tVolume (bbl)" << '\n';
    for (std::size_t i = 0; i < oil_tw.size(); ++i)
        std::cout << i << '\t' <<  oil_tw[i] << '\n';
    std::cout << "Oil TC: (qi = " << oil_tc.qi() / 365.25 << " bbl/d, Di = "
        << dca::convert_decline<dca::nominal, dca::secant_effective>(
                oil_tc.Di(), oil_tc.b()) * 100 << " sec. %/yr, b = "
        << oil_tc.b() << ")\n";
    std::cout << "Oil EUR: " << oil_eur / 1000 << " Mbbl\n";

    std::cout << "Gas Avg. Shift: " << gas_avg_shift << " months\n";
    std::cout << "Gas Type Well:\nMonth\tVolume (mcf)" << '\n';
    for (std::size_t i = 0; i < gas_tw.size(); ++i)
        std::cout << i << '\t' <<  gas_tw[i] << '\n';
    std::cout << "Gas TC: (qi = " << gas_tc.qi() / 365.25 << " mcf/d, Di = "
        << dca::convert_decline<dca::nominal, dca::secant_effective>(
                gas_tc.Di(), gas_tc.b()) * 100 << " sec. %/yr, b = "
        << gas_tc.b() << ")\n";
    std::cout << "Gas EUR: " << gas_eur / 1000 << " MMscf\n";

    std::cout << "6:1 BOE EUR: " << (oil_eur + gas_eur / 6) / 1000 << " Mboe\n";
}

dataset read_delimited(std::istream& is, char delim)
{
    dataset result;
    std::unordered_map<std::size_t, std::string> columns;
    std::string line;
    if (!std::getline(is, line))
        return result;

    std::size_t i = 0;
    for_delimited(line.begin(), line.end(), delim,
            [&](decltype(line.begin()) b, decltype(line.begin()) e) {
                columns[i++] = std::string(b, e);
            });

    while (std::getline(is, line)) {
        i = 0;
        for_delimited(line.begin(), line.end(), delim,
                [&](decltype(line.begin()) b, decltype(line.begin()) e) {
                    result[columns[i++]].emplace_back(b, e);
                });
        while (i < columns.size())
            result[columns[i++]].emplace_back();
    }

    return result;
}

template<class F>
F foreach_well(const dataset& data, F fn, std::string id_field)
{
    const auto& id = data.at(id_field);

    std::size_t begin_rec = 0, end_rec = 0;
    for (std::size_t i = 0; i < id.size(); ++i) {
        if (id[i] != id[begin_rec] || i == id.size() - 1) {
            if (i == id.size() - 1)
                end_rec = i;

            dataset well;
            std::for_each(data.begin(), data.end(),
                    [&](const std::pair<std::string,
                        std::vector<std::string>>& column)
                    {
                        well[column.first] = std::vector<std::string>(
                            column.second.data() + begin_rec,
                            column.second.data() + end_rec + 1);
                    }
            );
            fn(well);

            begin_rec = i;
        }
        end_rec = i;
    }

    return fn;
}
