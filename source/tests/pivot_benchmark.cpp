#include "tests.hpp"

#ifdef CMR_DO_TESTS

#include "lp_interface.hpp"
#include "core_lp.hpp"
#include "datagroups.hpp"
#include "util.hpp"

#include <catch.hpp>

#include <vector>
#include <string>
#include <iostream>
#include <iostream>

using std::min;
using std::max;
using std::abs;
using std::vector;

using std::string;
using std::cout;
using std::setprecision;

static double nd_time;
static double itlim_time;

static double nd_objval;
static double itlim_objval;

static int nd_itcount;
static int itlim_itcount;

SCENARIO ("Comparing pivot protocols",
          "[.LP][benchmark]") {
    vector<string> probs{
        "d2103",
        //"pr2392",
        // "pcb3038",
        // "rl5915",
        //"pla7397"
        // "usa13509"
    };


    for (string &prob : probs) {
        GIVEN ("The degree LP for " + prob) {
            CMR::Data::Instance inst("problems/" + prob + ".tsp", 99);
            CMR::Data::GraphGroup g_dat(inst);
            CMR::Data::BestGroup b_dat(inst, g_dat,
                "test_data/tours/" + prob + ".sol");
            CMR::LP::CoreLP core(g_dat, b_dat);

            double tourlen = b_dat.min_tour_value;


            THEN ("We can nd pivot to a new vector") {
                double t = CMR::util::zeit();
                core.nondegen_pivot(tourlen - CMR::Epsilon::Zero);
                nd_time = CMR::util::zeit() - t;
                nd_objval = core.get_objval();
                nd_itcount = core.it_count();
            }

            THEN ("We can single pivot to a new vector") {
                double t = CMR::util::zeit();
                int itcount = 0;
                while (core.get_objval() == tourlen) {
                    ++itcount;
                    core.single_pivot();
                }
                itlim_time = CMR::util::zeit() - t;
                itlim_objval = core.get_objval();
                itlim_itcount = itcount;

                cout << "\n\nInstance " << prob << "\n";
                cout << "piv times\t" << nd_time << "\t" << itlim_time << "\n"
                     << "piv vals\t" << nd_objval << "\t"
                     << itlim_objval << "\n"
                     << "piv counts\t" << nd_itcount << "\t" << itlim_itcount
                     << "\n";
            }
        }
    }
}

#endif
