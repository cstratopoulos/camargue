#include "config.hpp"

#ifdef CMR_DO_TESTS

#include "datagroups.hpp"
#include "karp.hpp"
#include "tooth.hpp"
#include "err_util.hpp"

#include <iostream>
#include <string>
#include <utility>
#include <vector>

extern "C" {
#include <concorde/INCLUDE/util.h>
}

#include <catch.hpp>

using std::cout;
using std::vector;
using std::unique_ptr;
using std::string;
using std::pair;


SCENARIO("Karp partition constructor tests",
	 "[.karp][valgrind]") {
  CMR::Graph::CoreGraph core_graph;
  CMR::Data::BestGroup b_dat;
  CMR::Data::SupportGroup s_dat;
  vector<double> lp_edges;
  CMR::Data::Instance inst;
  CMR::Data::KarpPartition kpart;

  REQUIRE_NOTHROW(CMR::Data::make_cut_test("problems/st70.tsp",
					    "test_data/tours/st70.sol",
					    "test_data/subtour_lp/st70.sub.x",
					    core_graph, b_dat, lp_edges, s_dat,
					    inst));
}

SCENARIO("Karp partitions for too small instances",
         "[.karp][tiny][.small]") {
  vector<string> probs{"blossom6", "fleisA9", "dantzig42"};
  for (string &fname : probs) {
    string
    probfile =
    ((fname == "dantzig42") ? "problems/" : "test_data/") + fname + ".tsp",
    solfile = "test_data/tours/" + fname + ".sol",
    subtourfile = "test_data/subtour_lp/" + fname + ".sub.x";
    CMR::Graph::CoreGraph core_graph;
    CMR::Data::BestGroup b_dat;
    CMR::Data::SupportGroup s_dat;
    vector<double> lp_edges;
    CMR::Data::Instance inst;
    CMR::Data::KarpPartition kpart;

    GIVEN("The TSPLIB instance " + fname) {
      THEN("The Karp partition is trivial bc there are too few nodes.") {
        REQUIRE_NOTHROW(CMR::Data::make_cut_test(probfile, solfile,
                                                 subtourfile, core_graph, b_dat,
                                                 lp_edges, s_dat, inst));
        int ncount = core_graph.node_count();

        REQUIRE_NOTHROW(kpart = CMR::Data::KarpPartition(inst));
        REQUIRE(kpart.num_parts() == 1);
      }
    }
  }
}

#endif //CMR_DO_TESTS
