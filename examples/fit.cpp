#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <cstddef>
#include <algorithm>
#include <cstdlib>

#include "dca/decline.hpp"
#include "dca/exponential.hpp"
#include "dca/hyperbolic.hpp"
#include "dca/hyptoexp.hpp"
#include "dca/bestfit.hpp"
#include "dca/production.hpp"

namespace params {
const std::string id_field = "Name";
static const std::string oil_field = "Oil";
static const std::string gas_field = "Gas";
static const auto d_final = dca::decline<dca::tangent_effective>(0.05);
static const double oil_el = 365.25;
static const double max_time = 30;
}

using dataset = std::unordered_map<std::string, std::vector<std::string>>;

dataset read_delimited(std::istream& is, char delim = '\t');

template<class F>
F foreach_well(const dataset& data, F fn,
        std::string id_field = params::id_field);

void process_well(const dataset& data);

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

    std::cout << params::id_field << '\t'
        << "OilEUR\tGasEUR\tBoeEUR\tOil.qi\tOil.Di\tOil.b\tOil.shift\t"
           "Gas.qi\tGas.Di\tGas.b\tGas.shift\n";
    foreach_well(data, process_well);
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

void process_well(const dataset& data)
{
    const auto& oil = data.at(params::oil_field);
    const auto& gas = data.at(params::gas_field);

    std::vector<double> oil_data(oil.size());
    std::transform(oil.begin(), oil.end(), oil_data.begin(),
            [](const std::string& d) {
                return std::strtod(d.c_str(), nullptr);
            });

    std::vector<double> gas_data(gas.size());
    std::transform(gas.begin(), gas.end(), gas_data.begin(),
            [](const std::string& d) {
                return std::strtod(d.c_str(), nullptr);
            });

    // strip leading zeros
    oil_data.erase(oil_data.begin(),
            std::find_if(oil_data.begin(), oil_data.end(),
                [](double p) { return p > 0.0; }));

    gas_data.erase(gas_data.begin(),
            std::find_if(gas_data.begin(), gas_data.end(),
                [](double p) { return p > 0.0; }));

    auto shifted_oil = dca::shift_to_peak(oil_data.begin(), oil_data.end());
    if (std::distance(std::get<0>(shifted_oil), oil_data.end()) < 3)
        return;
    auto oil_shift = std::distance(oil_data.begin(), std::get<0>(shifted_oil));
    auto oil_decline =
        dca::best_from_interval_volume<dca::arps_hyperbolic>(
            std::get<0>(shifted_oil), oil_data.end(), 0, 1.0 / 12.0);

    double t_eur;
    auto oil_eur = dca::eur(
            dca::arps_hyperbolic_to_exponential(
                oil_decline.qi(),
                oil_decline.Di(),
                oil_decline.b(),
                params::d_final
            ),
            params::oil_el, // bbl/yr = 1 bbl/day
            params::max_time, // years,
            &t_eur
    );

    auto shifted_gas = dca::shift_to_peak(gas_data.begin(), gas_data.end());
    auto gas_shift = std::distance(gas_data.begin(), std::get<0>(shifted_gas));
    auto gas_decline =
        dca::best_from_interval_volume<dca::arps_hyperbolic>(
            std::get<0>(shifted_gas), gas_data.end(), 0, 1.0 / 12.0);

    auto gas_eur = dca::arps_hyperbolic_to_exponential(
            gas_decline.qi(),
            gas_decline.Di(),
            gas_decline.b(),
            params::d_final
            ).cumulative(t_eur - (gas_shift - oil_shift));

    std::cout << data.at(params::id_field)[0] << '\t'
        << oil_eur / 1000 << '\t'
        << gas_eur / 1000 << '\t'
        << (oil_eur + gas_eur / 6) / 1000 << '\t'
        << oil_decline.qi() / 365.25 << '\t'
        << dca::convert_decline<dca::nominal, dca::secant_effective>(
                oil_decline.Di(), oil_decline.b()) << '\t'
        << oil_decline.b() << '\t'
        << oil_shift << '\t'
        << gas_decline.qi() / 365.25 << '\t'
        << dca::convert_decline<dca::nominal, dca::secant_effective>(
                gas_decline.Di(), gas_decline.b()) << '\t'
        << gas_decline.b() << '\t'
        << gas_shift << '\n';
}
