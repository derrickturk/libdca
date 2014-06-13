#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <cstddef>
#include <algorithm>
#include <cstdlib>

#include "exponential.hpp"
#include "hyperbolic.hpp"
#include "hyptoexp.hpp"
#include "bestfit.hpp"

using dataset = std::unordered_map<std::string, std::vector<std::string>>;

dataset read_delimited(std::istream& is, char delim = '\t');

template<class F>
F foreach_well(const dataset& data, F fn, std::string id_field = "PROPNUM");

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

    std::cout << "PROPNUM\tOilEUR\tGasEUR\tOil.qi\tOil.Di\tOil.b\t"
        "Gas.qi\tGas.Di\tGas.b\n";
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
    static const std::string id_field = "PROPNUM";
    static const std::string oil_field = "OIL";
    static const std::string gas_field = "GAS";

    const auto& oil = data.at(oil_field);
    const auto& gas = data.at(gas_field);

    std::vector<double> oil_data(oil.size());
    std::transform(oil.begin(), oil.end(), oil_data.begin(),
            [](const std::string& d) {
                return std::strtod(d.c_str(), nullptr);
            });

    auto peak_oil = std::max_element(oil_data.begin(), oil_data.end());
    if (std::distance(peak_oil, oil_data.end()) < 3)
        return;

    auto oil_decline =
        dca::best_from_interval_volume<dca::arps_hyperbolic>(
            peak_oil, oil_data.end(), 0, 1.0 / 12.0);

    double t_eur;
    auto oil_eur = dca::eur(
            dca::arps_hyperbolic_to_exponential(
                oil_decline.qi(),
                oil_decline.Di(),
                oil_decline.b(),
                dca::decline<dca::tangent_effective>(0.05)
            ),
            365.25, // bbl/yr = 1 bbl/day
            30, // years,
            &t_eur
    );

    std::vector<double> gas_data(gas.size());
    std::transform(gas.begin(), gas.end(), gas_data.begin(),
            [](const std::string& d) {
                return std::strtod(d.c_str(), nullptr);
            });

    auto gas_decline =
        dca::best_from_interval_volume<dca::arps_hyperbolic>(
            gas_data.begin() +
              std::distance(oil_data.begin(), peak_oil),
            gas_data.end(), 0, 1.0 / 12.0);

    auto gas_eur = dca::arps_hyperbolic_to_exponential(
            gas_decline.qi(),
            gas_decline.Di(),
            gas_decline.b(),
            dca::decline<dca::tangent_effective>(0.05)
            ).cumulative(t_eur);

    std::cout << data.at(id_field)[0] << '\t'
        << oil_eur << '\t'
        << gas_eur << '\t'
        << oil_decline.qi() / 365.25 << '\t'
        << dca::convert_decline<dca::nominal, dca::secant_effective>(
                oil_decline.Di(), oil_decline.b()) << '\t'
        << oil_decline.b() << '\t'
        << gas_decline.qi() / 365.25 << '\t'
        << dca::convert_decline<dca::nominal, dca::secant_effective>(
                gas_decline.Di(), gas_decline.b()) << '\t'
        << gas_decline.b() << '\n';
}
