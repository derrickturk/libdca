#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <iterator>
#include <unordered_map>
#include <cstddef>
#include <algorithm>
#include <cstdlib>
#include <functional>

#include "dca/exponential.hpp"
#include "dca/hyperbolic.hpp"
#include "dca/hyptoexp.hpp"
#include "dca/bestfit.hpp"
#include "dca/production.hpp"

namespace params {

const std::string id_field = "UID";
static const std::string oil_field = "Oil";
static const std::string gas_field = "Gas";
static const auto aggregation = dca::mean {};
static const auto d_final = dca::decline<dca::tangent_effective>(0.05);
static const double oil_el = 365.25;
static const double max_time = 30;

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
    using namespace std::placeholders;

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

    std::vector<std::pair<
        std::vector<double>::iterator,
        std::vector<double>::iterator>> oil_ranges;
    double oil_avg_shift = 0.0;
    for (auto& oil_rec : oil) {
        oil_rec.erase(std::find_if(oil_rec.rbegin(), oil_rec.rend(),
                    std::bind(std::not_equal_to<double> {}, 0.0, _1)).base(),
                oil_rec.end());
        oil_ranges.emplace_back(
                std::get<0>(dca::shift_to_peak(oil_rec.begin(), oil_rec.end())),
                oil_rec.end());
        oil_avg_shift += std::distance(oil_rec.begin(),
                oil_ranges.back().first);
    }
    oil_avg_shift /= oil_ranges.size();

    std::vector<std::pair<
        std::vector<double>::iterator,
        std::vector<double>::iterator>> gas_ranges;
    double gas_avg_shift = 0.0;
    for (auto& gas_rec : gas) {
        gas_rec.erase(std::find_if(gas_rec.rbegin(), gas_rec.rend(),
                    std::bind(std::not_equal_to<double> {}, 0.0, _1)).base(),
                gas_rec.end());
        gas_ranges.emplace_back(
                std::get<0>(dca::shift_to_peak(gas_rec.begin(), gas_rec.end())),
                gas_rec.end());
        gas_avg_shift += std::distance(gas_rec.begin(),
                gas_ranges.back().first);
    }
    gas_avg_shift /= gas_ranges.size();

    std::vector<double> oil_tw, gas_tw;
    dca::aggregate_production(oil_ranges.begin(), oil_ranges.end(),
            std::back_inserter(oil_tw), std::floor(oil_ranges.size() / 3),
            params::aggregation);
    dca::aggregate_production(gas_ranges.begin(), gas_ranges.end(),
            std::back_inserter(gas_tw), std::floor(gas_ranges.size() / 3),
            params::aggregation);

    auto oil_tc = dca::best_from_interval_volume<dca::arps_hyperbolic>(
            oil_tw.begin(), oil_tw.end(), 0, 1.0 / 12.0);
    auto gas_tc = dca::best_from_interval_volume<dca::arps_hyperbolic>(
            gas_tw.begin(), gas_tw.end(), 0, 1.0 / 12.0);

    std::vector<double> oil_forecast;
    dca::interval_volumes(oil_tc, std::back_inserter(oil_forecast),
            0.0, 1.0 / 12, oil_tw.size());

    std::vector<double> gas_forecast;
    dca::interval_volumes(gas_tc, std::back_inserter(gas_forecast),
            0.0, 1.0 / 12, gas_tw.size());

    double t_eur;
    double oil_eur = dca::eur(
            dca::arps_hyperbolic_to_exponential(
                oil_tc.qi(),
                oil_tc.Di(),
                oil_tc.b(),
                params::d_final
            ),
            params::oil_el, // bbl/yr = 1 bbl/day
            params::max_time, // years,
            &t_eur
    );

    double gas_eur = dca::arps_hyperbolic_to_exponential(
            gas_tc.qi(),
            gas_tc.Di(),
            gas_tc.b(),
            params::d_final
    ).cumulative(t_eur - (gas_avg_shift - oil_avg_shift));

    std::cout << "Oil Avg. Shift: " << oil_avg_shift << " months\n";
    std::cout << "Oil Type Well:\nMonth\tVolume (bbl)\tForecast (bbl)" << '\n';
    for (std::size_t i = 0; i < oil_tw.size(); ++i)
        std::cout << i << '\t' <<  oil_tw[i] << '\t' << oil_forecast[i] << '\n';
    std::cout << "Oil TC: (qi = " << oil_tc.qi() / 365.25 << " bbl/d, Di = "
        << dca::convert_decline<dca::nominal, dca::secant_effective>(
                oil_tc.Di(), oil_tc.b()) * 100 << " sec. %/yr, b = "
        << oil_tc.b() << ")\n";
    std::cout << "Oil EUR: " << oil_eur / 1000 << " Mbbl\n";

    std::cout << "Gas Avg. Shift: " << gas_avg_shift << " months\n";
    std::cout << "Gas Type Well:\nMonth\tVolume (mcf)\tForecast (mcf)" << '\n';
    for (std::size_t i = 0; i < gas_tw.size(); ++i)
        std::cout << i << '\t' <<  gas_tw[i] << '\t' << gas_forecast[i] << '\n';
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
