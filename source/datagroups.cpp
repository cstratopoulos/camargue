#include "datagroups.hpp"
#include "graph_io.hpp"
#include "util.hpp"
#include "err_util.hpp"

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
try : handle(util::make_unique<CCdatagroup>()), random_seed(seed) {
    if (CCutil_gettsplib(const_cast<char*>(fname.c_str()), &nodecount,
                         ptr()))
        throw runtime_error("CCutil_gettsplib failed.");
    
    cout << "Random seed " << seed << "\n";
    
    pname = fname.substr(fname.find_last_of("/") + 1);
    pname = pname.substr(0, pname.find_last_of("."));

    cout << std::fixed;

} catch (const exception &e) {
    cerr << e.what() << "\n";
    throw runtime_error("Instance constructor failed.");
 }

Instance::Instance(const int seed, const int ncount, const int gridsize)
try : handle(util::make_unique<CCdatagroup>()),
      nodecount(ncount),
      random_seed(seed),
      pname("r" + to_string(ncount) + "g" + to_string(gridsize)) {
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

  cout << std::fixed;

  
  cout << "Random problem, random seed " << seed << "\n";
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

    util::c_array_ptr edge_handle(elist);
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

    int kicks = 1000;
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
    cout << "LK run from best tour: " << min_tour_value << endl;

    for (int i = 0; i < perm.size(); ++i)
        perm[best_tour_nodes[i]] = i;

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

    for (int i = 0; i < best_tour_nodes.size(); ++i)
        perm[best_tour_nodes[i]] = i;

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

    cout << "Loaded and verified tour with length " << min_tour_value << endl;
} catch (const exception &e) {
    cerr << e.what() << "\n";
    throw runtime_error("BestGroup file constructor failed.");
}

void SupportGroup::reset(const int ncount, const vector<CMR::Edge> &edges,
                         const vector<double> &lp_x,
                         vector<int> &island)
{
    runtime_error err("Problem in SupportGroup::reset.");
    
    support_indices.clear();
    support_elist.clear();
    support_ecap.clear();

    integral = true;

    try {
        for (int i = 0; i < lp_x.size(); ++i)
            if (lp_x[i] >= CMR::Epsilon::Zero) {
                support_indices.push_back(i);
                support_ecap.push_back(lp_x[i]);
                support_elist.push_back(edges[i].end[0]);
                support_elist.push_back(edges[i].end[1]);

                if (lp_x[i] <= 1 - CMR::Epsilon::Zero)
                    integral = false;
            }
    } CMR_CATCH_PRINT_THROW("pushing back new support data", err);

    if (CMR::GraphUtils::build_s_graph(ncount, edges, support_indices, lp_x,
                                       &G_s))
        throw err;

    int icount = 0;
    connected = CMR::GraphUtils::connected(&G_s, &icount, island, 0);
}

bool SupportGroup::in_subtour_poly()
{
    if (!connected)
        return false;

    double cutval = 2;
    double rhs = 2.0 - CMR::Epsilon::Cut;

    if (CCcut_mincut(G_s.node_count, support_ecap.size(),
                     &support_elist[0], &support_ecap[0], &cutval,
                     NULL, NULL))
        throw runtime_error("CCcut_mincut failed in in_subtour_poly.");

    return cutval > rhs;
}

}
}
