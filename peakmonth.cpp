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

#include "exponential.hpp"
#include "hyperbolic.hpp"
#include "hyptoexp.hpp"
#include "bestfit.hpp"
#include "production.hpp"

namespace params {

const std::string id_field = "UID";
static const std::string oil_field = "Oil";
static const std::string gas_field = "Gas";
static const std::string month_field = "Month";
static const std::string api_field = "API";
static const std::string name_field = "Name";
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

    std::cout << "API\tName\tPeak Oil Month\t"
        "Peak Month Oil (bbl)\tPeak Oil Month Gas (mcf)\n";
    foreach_well(data, [&](const dataset& well) {
        auto&& oil_text = well.at(params::oil_field);
        std::vector<double> oil (oil_text.size());
        std::transform(oil_text.begin(), oil_text.end(), oil.begin(),
            [](const std::string& s) { return std::strtod(s.c_str(), nullptr); }
        );

        auto&& gas_text = well.at(params::gas_field);
        std::vector<double> gas (gas_text.size());
        std::transform(gas_text.begin(), gas_text.end(), gas.begin(),
            [](const std::string& s) { return std::strtod(s.c_str(), nullptr); }
        );

        auto peak = dca::shift_to_peak(oil.begin(), oil.end(), gas.begin());

        std::cout << well.at(params::api_field)[0] << '\t'
          << well.at(params::name_field)[0] << '\t'
          << well.at(params::month_field)[
               std::distance(oil.begin(), std::get<0>(peak))]
          << '\t' << *std::get<0>(peak) << '\t'
          << *std::get<1>(peak) << '\n';
    });
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
