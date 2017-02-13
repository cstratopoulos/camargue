#include "config.hpp"

#ifdef CMR_DO_TESTS

#include "lp_interface.hpp"
#include "core_lp.hpp"
#include "datagroups.hpp"
#include "util.hpp"

#include <catch.hpp>

#include <algorithm>
#include <array>
#include <vector>
#include <string>
#include <iostream>
#include <iostream>

using std::array;
using std::min;
using std::max;
using std::abs;
using std::vector;

using std::string;
using std::cout;
using std::endl;
using std::setprecision;

template <typename T>
using Triple = std::array<T, 3>;
using RepTuple = std::tuple<string, Triple<int>, Triple<double>, int>;

static vector<RepTuple> table_entries{
    RepTuple("d2103", {{0,0,0}}, {{0.0, 0.0, 0.0}}, 2103),
    RepTuple("fl3795", {{0,0,0}}, {{0.0, 0.0, 0.0}}, 3795),
    // RepTuple("fnl4461", {{0,0,0}}, {{0.0, 0.0, 0.0}}, 4461),
    // RepTuple("pcb3038", {{0,0,0}}, {{0.0, 0.0, 0.0}}, 3038),
    // RepTuple("pla7397", {{0,0,0}}, {{0.0, 0.0, 0.0}}, 7397),
    // RepTuple("pr2392", {{0,0,0}}, {{0.0, 0.0, 0.0}}, 2392),
    // RepTuple("rl5915", {{0,0,0}}, {{0.0, 0.0, 0.0}}, 5915),
    // RepTuple("rl5934", {{0,0,0}}, {{0.0, 0.0, 0.0}}, 5934),
    // RepTuple("u2152", {{0,0,0}}, {{0.0, 0.0, 0.0}}, 2152),
    // RepTuple("u2319", {{0,0,0}}, {{0.0, 0.0, 0.0}}, 2319),
    // RepTuple("brd14051", {{0,0,0}}, {{0.0, 0.0, 0.0}}, 14051),
    // RepTuple("rl11849", {{0,0,0}}, {{0.0, 0.0, 0.0}}, 11849),
    };

SCENARIO ("Comparing pivot protocols as optimizers",
          "[LP][Relaxation][single_pivot][nondegen_pivot][benchmark]") {
    using namespace CMR;
    namespace Eps = Epsilon;

    for (RepTuple &te : table_entries) {
        string prob = std::get<0>(te);
        Triple<int> &piv_counts = std::get<1>(te);
        Triple<double> &piv_times = std::get<2>(te);
        int ncount = std::get<3>(te);
        GIVEN ("The degree LP for " + prob) {
            Data::Instance inst("problems/" + prob + ".tsp", 99);
            Data::GraphGroup g_dat(inst);
            Data::BestGroup b_dat(inst, g_dat,
                                  "test_data/tours/" + prob + ".sol");
            LP::CoreLP core(g_dat, b_dat);

            double tourlen = b_dat.min_tour_value;
            REQUIRE(core.get_objval() == tourlen);

            THEN ("We can primal opt the degree LP") {
                cout << "Testing primal opt..." << endl;
                double t = util::zeit();
                core.primal_opt();
                piv_times[0] = util::zeit() - t;
                piv_counts[0] = core.it_count();                
            }

            THEN ("We can primal opt one nd pivot at a time") {
                cout << "Testing nd opt..." << endl;
                double t = util::zeit();
                int nd_itcount = 0;
                while (!core.dual_feas()) {
                    double objval = core.get_objval();
                    core.nondegen_pivot(objval - Eps::Zero);
                    nd_itcount += core.it_count();
                }
                piv_times[1] = util::zeit() - t;
                piv_counts[1] = nd_itcount;
            }

            THEN ("We can primal opt a single pivot at a time") {
                cout << "Testing itlim opt..." << endl;
                double t = util::zeit();
                int it_itcount = 0;
                while (!core.dual_feas()) {
                    ++it_itcount;
                    core.one_primal_pivot();
                }
                piv_times[2] = util::zeit() - t;
                piv_counts[2] = it_itcount;
            }
        }
    }
    cout << "\n";

    THEN ("Report the results") {
        std::sort(table_entries.begin(), table_entries.end(),
                  [](RepTuple r, RepTuple t)
                  { return std::get<3>(r) < std::get<3>(t); });

        cout << "Instance\tCPU Time (scaled)\tIteration count\n";
        for (RepTuple &te : table_entries) {
            string prob = std::get<0>(te);
            Triple<int> &piv_counts = std::get<1>(te);
            Triple<double> &piv_times = std::get<2>(te);
            using DescPair = std::pair<string, int>;
            std::vector<DescPair>
            reports{DescPair("Primal opt", 0),
                DescPair("Nondeg opt", 1),
                DescPair("Itlim opt", 2)};
            cout << prob << "\n";
            for (auto dp : reports) {
                string protocol = dp.first;
                int index = dp.second;
                printf("%s\t%.2f\t%d\n",
                       protocol.c_str(),
                       (piv_times[index] / piv_times[0]),
                       piv_counts[index]);
            }
        }
    }
    
}

#endif
