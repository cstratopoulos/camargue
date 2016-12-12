#include "datagroups.hpp"
#include "graph_io.hpp"
#include "util.hpp"

#include <algorithm>
#include <unordered_map>
#include <iostream>
#include <iomanip>
#include <utility>
#include <vector>
#include <stdexcept>

#include <cmath>

extern "C" {
#include <concorde/INCLUDE/linkern.h>
#include <concorde/INCLUDE/edgegen.h>
}

using std::cout;
using std::cerr;
using std::endl;
using std::setprecision;

using std::string;
using std::to_string;

using std::vector;
using std::unique_ptr;

using std::min;
using std::max;

using std::exception;
using std::runtime_error;
using std::logic_error;

namespace CMR {
namespace Data {

Instance::Instance() noexcept : handle(nullptr) {}

Instance::Instance(const string &fname, const int seed)
try : handle(CMR::make_unique<CCdatagroup>()),
      random_seed(seed) {
    
    if (CCutil_gettsplib(const_cast<char*>(fname.c_str()), &nodecount,
                         ptr()))
        throw runtime_error("CCutil_gettsplib failed.");
    
    pname = fname.substr(fname.find_last_of("/") + 1);
    pname = pname.substr(0, pname.find_last_of("."));
} catch (const exception &e) {
    cerr << e.what() << "\n";
    throw runtime_error("Instance constructor failed.");
 }

Instance::Instance(const int seed, const int ncount, const int gridsize)
try : handle(CMR::make_unique<CCdatagroup>()),
      nodecount(ncount),
      random_seed(seed),
      pname("rand_s" + to_string(seed) + "n" + to_string(ncount) + "g"
            + to_string(gridsize)) {
  if (ncount <= 0)
      throw logic_error("Specified bad ncount.");
  
  if (gridsize <= 0)
      throw logic_error("Specified bad gridsize.");
    
  CCrandstate rstate;
  int allow_dups = 1;
  int binary_in = 0;

  int tmp_ncount = ncount;
  int tmp_gridsize = gridsize;

  CCutil_sprand(seed, &rstate);
  if (CCutil_getdata((char *) NULL, binary_in, CC_EUCLIDEAN,
		    &tmp_ncount, ptr(), tmp_gridsize, allow_dups, &rstate))
    throw runtime_error("CCutil_getdata failed.");
} catch (const exception &e) {
  cerr << e.what() << "\n";
  throw runtime_error("Instance constructor failed.");
 }

Instance::Instance(Instance &&I) noexcept :
    handle(std::move(I.handle)), nodecount(I.nodecount),
    random_seed(I.random_seed), pname(I.pname)
{}

Instance &Instance::operator=(Instance &&I) noexcept {
  handle = std::move(I.handle);
  nodecount = I.nodecount;
  random_seed = I.random_seed;
  pname = I.pname;
  return *this;
}

GraphGroup::GraphGroup(const Instance &inst)
try     
{
    int ncount = inst.node_count();
    m_graph.node_count = ncount;

    CCedgegengroup plan;
    CCrandstate rstate;

    CCutil_sprand(inst.seed(), &rstate);
    CCedgegen_init_edgegengroup(&plan);

    plan.linkern.count = 9;
    plan.linkern.quadnearest = 2;
    plan.linkern.greedy_start = 0;
    plan.linkern.nkicks = (ncount / 100) + 1;

    int *elist = (int *) NULL;

    if(CCedgegen_edges(&plan, ncount, inst.ptr(), NULL,
                       &(m_graph.edge_count), &elist, 1, &rstate))
        throw runtime_error("CCedgegen_edges failed.");

    c_array_ptr edge_handle(elist);
    int ecount = m_graph.edge_count;

    for (int i = 0; i < ecount; ++i) {
        int end0 = min(elist[2 * i], elist[(2 * i) + 1]);
        int end1 = max(elist[2 * i], elist[(2 * i) + 1]);

        m_graph.edges.emplace_back(end0, end1,
                                   CCutil_dat_edgelen(end0, end1, inst.ptr()));
        m_graph.edge_lookup.emplace(IntPair(end0, end1), i);
    }

    island.resize(ncount);
    delta.resize(ecount, 0);
    edge_marks.resize(ncount, 0);
} catch (const exception &e) {
    cerr << e.what() << "\n";
    throw runtime_error("GraphGroup constructor failed.");
}

BestGroup::BestGroup(const Instance &inst, GraphGroup &graph_data) try :
    best_tour_edges(std::vector<int>(graph_data.m_graph.edge_count, 0)),
    best_tour_nodes(std::vector<int>(graph_data.m_graph.node_count)),
    perm(best_tour_nodes.size()),
    min_tour_value(DoubleMax)
{

    CCrandstate rstate;
    CCutil_sprand(inst.seed(), &rstate);

    Graph &graph = graph_data.m_graph;
    int ncount = graph.node_count;
    int ecount = graph.edge_count;
    vector<int> elist;

    for(Edge &e : graph.edges) {
        elist.push_back(e.end[0]);
        elist.push_back(e.end[1]);
    }

    int kicks = 500;
    int trials = 10;
    
    if (CClinkern_tour(ncount, inst.ptr(), ecount, &elist[0], ncount, kicks,
                       (int *) NULL, &best_tour_nodes[0], &min_tour_value,
                       1, 0.0, 0.0, (char *) NULL, CC_LK_GEOMETRIC_KICK,
                       &rstate))
        throw runtime_error("CClinkern_tour failed.");

    cout << "LK initial run: " << min_tour_value << ". Performing "
         << trials << " more trials. (";

    vector<int> cyc(ncount);
    double tourlen(DoubleMax);

    for (int i = 0; i < trials; ++i) {
        if (CClinkern_tour(ncount, inst.ptr(), ecount, &elist[0], ncount,
                           kicks, (int *) NULL, &cyc[0], &tourlen, 1, 0.0, 0.0,
                           (char *) NULL, CC_LK_GEOMETRIC_KICK, &rstate))
            throw runtime_error("CClinkern_tour failed.");

        if (tourlen < min_tour_value) {
            best_tour_nodes = cyc;
            min_tour_value = tourlen;
            cout << "!";
        } else
            cout << ".";
        cout.flush();
    }

    cout << ")\n";

    if (CClinkern_tour(ncount, inst.ptr(), ecount, &elist[0], ncount,
                       2 * kicks, &best_tour_nodes[0], &cyc[0], &min_tour_value,
                       1, 0.0, 0.0, (char *) NULL, CC_LK_GEOMETRIC_KICK,
                       &rstate))
        throw runtime_error("CClinkern_tour failed.");

    best_tour_nodes = cyc;
    cout << "LK run from best tour: " << min_tour_value << "\n";

    for (int i : best_tour_nodes)
        perm[i] = i;

    vector<int> &delta = graph_data.delta;

    for (int i = 0; i < ncount; ++i) {
        int end0 = min(best_tour_nodes[i], best_tour_nodes[(i + 1) % ncount]);
        int end1 = max(best_tour_nodes[i], best_tour_nodes[(i + 1) % ncount]);
        IntPair find_pair(end0, end1);
        IntPairMap::iterator edge_it = graph.edge_lookup.find(find_pair);

        if (edge_it == graph.edge_lookup.end()) {
            Edge e(end0, end1, CCutil_dat_edgelen(end0, end1, inst.ptr()));

            graph.edges.push_back(e);
            graph.edge_lookup[find_pair] = graph.edges.size() - 1;
            graph.edge_count += 1;
            best_tour_edges.push_back(0);
            delta.push_back(0);
            edge_it = graph.edge_lookup.find(find_pair);
        }

        int edge_index = edge_it->second;
        best_tour_edges[edge_index] = 1;
    }

    if ((ncount % 2) == 0) {
        int end0 = fmin(best_tour_nodes[0], best_tour_nodes[ncount - 2]);
        int end1 = fmax(best_tour_nodes[0], best_tour_nodes[ncount - 2]);
        IntPair find_pair(end0, end1);
        IntPairMap::iterator edge_it = graph.edge_lookup.find(find_pair);

        if (edge_it == graph.edge_lookup.end()) {
            Edge e(end0, end1, CCutil_dat_edgelen(end0, end1, inst.ptr()));

            graph.edges.push_back(e);
            graph.edge_lookup[find_pair] = graph.edges.size() - 1;
            graph.edge_count += 1;
            best_tour_edges.push_back(0);
            delta.push_back(0);
        }
    }
} catch (const exception &e) {
    cerr << e.what() << "\n";
    throw runtime_error("BestGroup LK constructor failed.");
}

BestGroup::BestGroup(const Instance &inst, GraphGroup &graph_data,
          const std::string &tourfile) try :
    best_tour_edges(std::vector<int>(graph_data.m_graph.edge_count, 0)),
    best_tour_nodes(std::vector<int>(graph_data.m_graph.node_count)),
    perm(best_tour_nodes.size()),
    min_tour_value(DoubleMax)
{
    Graph &graph = graph_data.m_graph;
    int ncount = graph.node_count;

    get_tour_nodes(ncount, best_tour_nodes, tourfile);

    for (int i : best_tour_nodes)
        perm[i] = i;

    vector<int> &delta = graph_data.delta;

    for (int i = 0; i < ncount; ++i) {
        int end0 = min(best_tour_nodes[i], best_tour_nodes[(i + 1) % ncount]);
        int end1 = max(best_tour_nodes[i], best_tour_nodes[(i + 1) % ncount]);
        IntPair find_pair(end0, end1);
        IntPairMap::iterator edge_it = graph.edge_lookup.find(find_pair);

        if (edge_it == graph.edge_lookup.end()) {
            Edge e(end0, end1, CCutil_dat_edgelen(end0, end1, inst.ptr()));

            graph.edges.push_back(e);
            graph.edge_lookup[find_pair] = graph.edges.size() - 1;
            graph.edge_count += 1;
            best_tour_edges.push_back(0);
            delta.push_back(0);
            edge_it = graph.edge_lookup.find(find_pair);
        }

        int edge_index = edge_it->second;
        best_tour_edges[edge_index] = 1;
    }

    if ((ncount % 2) == 0) {
        int end0 = fmin(best_tour_nodes[0], best_tour_nodes[ncount - 2]);
        int end1 = fmax(best_tour_nodes[0], best_tour_nodes[ncount - 2]);
        IntPair find_pair(end0, end1);
        IntPairMap::iterator edge_it = graph.edge_lookup.find(find_pair);

        if (edge_it == graph.edge_lookup.end()) {
            Edge e(end0, end1, CCutil_dat_edgelen(end0, end1, inst.ptr()));

            graph.edges.push_back(e);
            graph.edge_lookup[find_pair] = graph.edges.size() - 1;
            graph.edge_count += 1;
            best_tour_edges.push_back(0);
            delta.push_back(0);
        }
    }

    min_tour_value = 0;

    for (int i = 0; i < best_tour_edges.size(); ++i)
        if (best_tour_edges[i] == 1)
            min_tour_value += graph.edges[i].len;
} catch (const exception &e) {
    cerr << e.what() << "\n";
    throw runtime_error("BestGroup file constructor failed.");
}


LPGroup::LPGroup(const Graph &m_graph, CMR::LP::Prefs &_prefs,
			   const vector<int> &perm) {
  int rval = 0;
  int cmatbeg = 0, num_vars = 1, num_non_zero = 2;
  double coefficients[2] = {1.0, 1.0};
  double lower_bound = 0.0;
  double upper_bound = 1.0;

  
  //Build the basic LP
  rval = CMRlp_init (&m_lp);
  if (rval) goto CLEANUP;

  rval = CMRlp_create (&m_lp, "subtour");
  if (rval) goto CLEANUP;
	  

  /* Build a row for each degree equation */
  for (int i = 0; i < m_graph.node_count; i++) {
    rval = CMRlp_new_row (&m_lp, 'E', 2.0);
    if (rval) goto CLEANUP;
  }

  /* Build a column for each edge of the Graph */
  for (int j = 0; j < m_graph.edge_count; j++) {
    int *nodes = (int*)m_graph.edges[j].end;
    double objective_val = (double)m_graph.edges[j].len;
    rval = CMRlp_addcols (&m_lp, num_vars, num_non_zero, &objective_val,
			   &cmatbeg, nodes, coefficients, &lower_bound,
			   &upper_bound);
    if (rval) goto CLEANUP;
  }

  prefs = _prefs;

  try {
    m_lp_edges.resize(m_graph.edge_count);
    old_colstat.resize(m_graph.edge_count, CPX_AT_LOWER);
    old_rowstat.resize(m_graph.node_count, CPX_AT_LOWER);
    frac_colstat.resize(m_graph.edge_count);
    frac_rowstat.resize(m_graph.edge_count);
  } catch (const std::bad_alloc &) {
    rval = 1; CMR_GOTO_CLEANUP("Problem allocating LP vectors, ");
  }

 CLEANUP:
  if (rval) {
    cerr << "LPGroup constructor failed\n";
    throw 1;
  }
}

}
}
