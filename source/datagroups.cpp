#include "datagroups.hpp"
#include "io_util.hpp"
#include "util.hpp"
#include "err_util.hpp"

#include <algorithm>

#include <iostream>
#include <iomanip>
#include <limits>
#include <utility>
#include <vector>
#include <stdexcept>

#include <cmath>

extern "C" {
#include <concorde/INCLUDE/cut.h>
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
using std::logic_error;
using std::runtime_error;

namespace CMR {

namespace Eps = Epsilon;
constexpr double DoubleMax = std::numeric_limits<double>::max();

namespace Data {

Instance::Instance() noexcept { CCutil_init_datagroup(&dat); }


/**
 * @param[in] fname a path to a TSPLIB file specified from the executable
 * directory.
 * @param[in] seed the random seed to be used throughout.
 */
Instance::Instance(const string &fname, const int seed)
try : random_seed(seed) {
    CCutil_init_datagroup(&dat);

    if (CCutil_gettsplib(const_cast<char*>(fname.c_str()), &nodecount,
                         ptr()))
        throw runtime_error("CCutil_gettsplib failed.");

    cout << "Random seed " << seed << endl;

    pname = fname.substr(fname.find_last_of("/") + 1);
    pname = pname.substr(0, pname.find_last_of("."));

    cout << std::fixed;

} catch (const exception &e) {
    cerr << e.what() << "\n";
    throw runtime_error("Instance constructor failed.");
}

/**
 * @param[in] seed the random seed used to generate the problem, and for all
 * other later random computations.
 * @param[in] ncount the number of nodes in the problem
 * @param[in] gridsize the nodes will be generated from a square grid with this
 * side length.
 */
Instance::Instance(const int seed, const int ncount, const int gridsize)
try : nodecount(ncount), random_seed(seed),
      pname("r" + to_string(ncount) + "g" + to_string(gridsize)
            + "s" + to_string(seed)) {
  if (ncount <= 2)
      throw logic_error("Specified bad ncount.");

  if (gridsize <= 0)
      throw logic_error("Specified bad gridsize.");

  CCrandstate rstate;
  int allow_dups = 1;
  int binary_in = 0;

  int tmp_ncount = ncount;
  int tmp_gridsize = gridsize;

  CCutil_sprand(seed, &rstate);

  CCutil_init_datagroup(&dat);

  if (CCutil_getdata((char *) NULL, binary_in, CC_EUCLIDEAN,
		    &tmp_ncount, ptr(), tmp_gridsize, allow_dups, &rstate))
    throw runtime_error("CCutil_getdata failed.");

  cout << std::fixed;


  cout << "Random problem, random seed " << seed << endl;
} catch (const exception &e) {
  cerr << e.what() << "\n";
  throw runtime_error("Instance constructor failed.");
}

/**
 * This constructor generates a sparse Instance which contains precisely the
 * edges and lengths specified by \p elist and \p elen.
 * @param probname the name of the problem.
 * @param seed the random seed to store.
 * @param ncount the number of nodes.
 * @param elist a node-node elist of the sparse edges.
 * @param elen the edge capacities on the sparse edges.
 * @param default_len the length to be assigned to edges not in the graph.
 * If a non-positive value is selected, a large length on the order of
 * \p ncount times the largest length in \p elen will be chosen.
 */
Instance::Instance(const string &probname, int seed, int ncount,
                   vector<int> &elist, vector<int> &elen, int default_len)
try
    : nodecount(ncount), random_seed(seed), pname(probname)
{
    if (default_len <= 0)
        cout << "\tCreating sparse Instance...";
    if (CCutil_graph2dat_sparse(ncount, elen.size(), &elist[0], &elen[0],
                                default_len, ptr()))
        throw runtime_error("CCutil_getdata failed.");

    cout << std::fixed;
} catch (const exception &e) {
    cerr << e.what() << "\n";
    throw runtime_error("Instance constructor failed.");
}

Instance::Instance(Instance &&I) noexcept :
    nodecount(I.nodecount), random_seed(I.random_seed), pname(I.pname)
{
    CCutil_freedatagroup(&dat);
    dat = I.dat;

    CCutil_init_datagroup(&I.dat);
    I.nodecount = 0;
    I.random_seed = 0;
    I.pname.clear();
}

Instance& Instance::operator=(Instance &&I) noexcept
{
  nodecount = I.nodecount;
  random_seed = I.random_seed;
  pname = I.pname;

  CCutil_freedatagroup(&dat);
  dat = I.dat;

  CCutil_init_datagroup(&I.dat);
  I.nodecount = 0;
  I.random_seed = 0;
  I.pname.clear();

  return *this;
}

Instance::~Instance() { CCutil_freedatagroup(&dat); }

double Instance::tour_length(const vector<int> &tour_nodes) const
{
    if (tour_nodes.size() != nodecount) {
        cerr << "Tour vector size " << tour_nodes.size() << " with nodecount "
             << nodecount << endl;
        throw runtime_error("Size mismatch in Instance::tour_length");
    }

    double result = 0.0;

    for (int i = 0; i < nodecount; ++i)
        result += edgelen(tour_nodes[i], tour_nodes[(i + 1) % nodecount]);

    return result;
}

}

namespace Graph {

CoreGraph::CoreGraph(const Data::Instance &inst, Graph::EdgePlan edge_plan) try
    : nodecount(inst.node_count())
{
    int ncount = nodecount;
    int ecount = 0;

    CCedgegengroup plan;
    CCrandstate rstate;
    int norm = inst.ptr()->norm;

    if (edge_plan == EdgePlan::Delaunay)
        if (norm != CC_EUCLIDEAN && norm != CC_EUCLIDEAN_CEIL) {
            cout << "Requested Delaunay for incompatible norm, overriding"
                 << endl;
            edge_plan = EdgePlan::Linkern;
        }

    CCutil_sprand(inst.seed(), &rstate);
    CCedgegen_init_edgegengroup(&plan);

    switch (edge_plan) {
    case EdgePlan::Linkern:
        cout << "Linkern edges\n";
        plan.linkern.count = 9;
        plan.linkern.quadnearest = 2;
        plan.linkern.greedy_start = 0;
        plan.linkern.nkicks = (ncount / 100) + 1;
        if (ncount < 100)
            plan.quadnearest = 2;
        break;
    case EdgePlan::Delaunay:
        cout << "Delaunay edges\n";
        plan.delaunay = 1;
        break;
    default:
        throw logic_error("Unimplemented EdgePlan case.");
    }

    int *elist = (int *) NULL;

    if(CCedgegen_edges(&plan, ncount, inst.ptr(), NULL,
                       &ecount, &elist, 1, &rstate))
        throw runtime_error("CCedgegen_edges failed.");

    util::c_array_ptr<int> edge_handle(elist);

    edges.reserve(ecount);

    for (int i = 0; i < ecount; ++i) {
        int e0 = elist[2 * i];
        int e1 = elist[(2 * i) + 1];
        edges.emplace_back(e0, e1, inst.edgelen(e0, e1));
    }

    adj_list = AdjList(ncount, edges);

    cout << "Initialized with " << ecount << " edges." << endl;
} catch (const exception &e) {
    cerr << e.what() << "\n";
    throw runtime_error("CoreGraph Instance constructor failed.");
}

/**
 * @param ncount the number of nodes in the graph.
 * @param ecount the number of edges.
 * @param elist a node-node elist of length \p ecount.
 * @param edgelen a function determining the length of edges in \p elist.
 */
CoreGraph::CoreGraph(int ncount, int ecount, const int *elist,
                     const std::function<double(int, int)> edgelen) try
    : nodecount(ncount)
{
    edges.reserve(ecount);

    for (int i = 0; i < ecount; ++i) {
        int e0 = elist[2 * i];
        int e1 = elist[(2 * i) + 1];
        edges.emplace_back(e0, e1, edgelen(e0, e1));
    }

    adj_list = AdjList(ncount, edges);
} catch (const exception &e) {
    cerr << e.what() << "\n";
    throw runtime_error("CoreGraph constructor failed.");
}

/**
 * @param tour_nodes the tour sequence to represent.
 * @param edgelen the edge length function to determine edge weights.
 */
CoreGraph::CoreGraph(const vector<int> &tour_nodes,
                     const std::function<double(int, int)> edgelen) try
    : nodecount(tour_nodes.size())
{
    edges.reserve(nodecount);

    for (int i = 0; i < nodecount; ++i) {
        int e0 = tour_nodes[i];
        int e1 = tour_nodes[(i + 1) % nodecount];

        edges.emplace_back(e0, e1, edgelen(e0, e1));
    }

    adj_list = AdjList(nodecount, edges);
} catch (const exception &e) {
    cerr << e.what() << "\n";
    throw runtime_error("CoreGraph constructor failed.");
}

int CoreGraph::find_edge_ind(int end0, int end1) const
{
    const AdjObj *adj_ptr = adj_list.find_edge(end0, end1);
    if (adj_ptr == nullptr)
        return -1;
    int result = adj_ptr->edge_index;

    if (result > edge_count()) {
        cerr << end0 << ", " << end1 << " found with index "
             << result << " but " << edge_count() << " edges" << endl;
        throw runtime_error("Result out of bounds in CoreGraph::find_edge_ind");
    }

    return adj_ptr->edge_index;
}

void CoreGraph::add_edge( int end0, int end1, int len )
{
    if (find_edge_ind(end0, end1) != -1)
        return;

    int new_ind = edge_count();

    edges.emplace_back(end0, end1, len);
    adj_list.add_edge(end0, end1, new_ind, len);
}

void CoreGraph::add_edge(const Edge &e)
{
    if (find_edge_ind(e.end[0], e.end[1]) != -1)
        return;

    int new_ind = edge_count();
    edges.push_back(e);
    adj_list.add_edge(e.end[0], e.end[1], new_ind, e.len);
}

void CoreGraph::remove_edges()
{
    edges.erase(std::remove_if(edges.begin(), edges.end(),
                               [](const Edge &e) { return e.removable; }),
                edges.end());

    adj_list = AdjList(node_count(), edges);
}

}

namespace Data {

BestGroup::BestGroup(const Instance &inst, Graph::CoreGraph &core_graph) try :
    best_tour_edges(std::vector<int>(core_graph.edge_count(), 0)),
    best_tour_nodes(std::vector<int>(core_graph.node_count())),
    perm(best_tour_nodes.size()),
    min_tour_value(DoubleMax)
{

    CCrandstate rstate;
    CCutil_sprand(inst.seed(), &rstate);

    int ncount = core_graph.node_count();
    int ecount = core_graph.edge_count();
    vector<int> elist;

    for(const Graph::Edge &e : core_graph.get_edges()) {
        elist.push_back(e.end[0]);
        elist.push_back(e.end[1]);
    }

    int kicks = 1000;
    int trials = 10;
    int kicktype = ((ncount < 10000) ? CC_LK_RANDOM_KICK :
                    CC_LK_GEOMETRIC_KICK);

    if (CClinkern_tour(ncount, inst.ptr(), ecount, &elist[0], ncount, kicks,
                       (int *) NULL, &best_tour_nodes[0], &min_tour_value,
                       1, 0.0, 0.0, (char *) NULL, kicktype,
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

    cout << ")" << endl;

    if (CClinkern_tour(ncount, inst.ptr(), ecount, &elist[0], ncount,
                       2 * kicks, &best_tour_nodes[0], &cyc[0], &min_tour_value,
                       1, 0.0, 0.0, (char *) NULL, CC_LK_GEOMETRIC_KICK,
                       &rstate))
        throw runtime_error("CClinkern_tour failed.");

    best_tour_nodes = cyc;
    cout << "LK run from best tour: " << min_tour_value << endl;

    for (int i = 0; i < perm.size(); ++i)
        perm[best_tour_nodes[i]] = i;

    for (int i = 0; i < ncount; ++i) {
        int e0 = best_tour_nodes[i];
        int e1 = best_tour_nodes[(i + 1) % ncount];

        int find_ind = core_graph.find_edge_ind(e0, e1);

        if (find_ind == -1) {
            core_graph.add_edge(e0, e1, inst.edgelen(e0, e1));
            find_ind = core_graph.find_edge_ind(e0, e1);

            best_tour_edges.push_back(0);
        }

        best_tour_edges[find_ind] = 1;
    }

    if ((ncount % 2) == 0) {
        int e0 = best_tour_nodes[0];
        int e1 = best_tour_nodes[ncount - 2];

        int find_ind = core_graph.find_edge_ind(e0, e1);

        if (find_ind == -1) {
            core_graph.add_edge(e0, e1, inst.edgelen(e0, e1));
            find_ind = core_graph.find_edge_ind(e0, e1);

            best_tour_edges.push_back(0);
        }
    }
} catch (const exception &e) {
    cerr << e.what() << "\n";
    throw runtime_error("BestGroup LK constructor failed.");
}

BestGroup::BestGroup(const Instance &inst, Graph::CoreGraph &core_graph,
                     const std::string &tourfile) try :
    best_tour_edges(std::vector<int>(core_graph.edge_count(), 0)),
    best_tour_nodes(std::vector<int>(core_graph.node_count())),
    perm(best_tour_nodes.size()),
    min_tour_value(DoubleMax)
{
    int ncount = core_graph.node_count();

    util::get_tour_nodes(ncount, best_tour_nodes, tourfile);

    for (int i = 0; i < best_tour_nodes.size(); ++i)
        perm[best_tour_nodes[i]] = i;

    for (int i = 0; i < ncount; ++i) {
        int e0 = best_tour_nodes[i];
        int e1 = best_tour_nodes[(i + 1) % ncount];

        int find_ind = core_graph.find_edge_ind(e0, e1);

        if (find_ind == -1) {
            core_graph.add_edge(e0, e1, inst.edgelen(e0, e1));
            find_ind = core_graph.find_edge_ind(e0, e1);

            best_tour_edges.push_back(0);
        }

        best_tour_edges[find_ind] = 1;
    }


    if ((ncount % 2) == 0) {
        int e0 = best_tour_nodes[0];
        int e1 = best_tour_nodes[ncount - 2];

        int find_ind = core_graph.find_edge_ind(e0, e1);

        if (find_ind == -1) {
            core_graph.add_edge(e0, e1, inst.edgelen(e0, e1));
            find_ind = core_graph.find_edge_ind(e0, e1);

            best_tour_edges.push_back(0);
        }
    }

    min_tour_value = 0;

    for (int i = 0; i < best_tour_edges.size(); ++i)
        if (best_tour_edges[i] == 1)
            min_tour_value += core_graph.get_edge(i).len;

    cout << "Loaded and verified tour with length " << min_tour_value << endl;
} catch (const exception &e) {
    cerr << e.what() << "\n";
    throw runtime_error("BestGroup file constructor failed.");
}

/**
 * @param[in] edges the edges of the CoreGraph relative to which this
 * SupportGroup is being constructed.
 * @param[in] lp_x the lp solution vector.
 * @param[out] island a list of nodes in a connected component found in
 * the dfs.
 * @param[in] ncount the number of nodes in the graph.
 * Constructs a SupportGroup with respect to an lp solution vector, with edge
 * capacities corresponding to the lp values in \p lp_x. Computes basic
 * connectivity info.
 */
SupportGroup::SupportGroup(const vector<Graph::Edge> &edges,
                           const vector<double> &lp_x,
                           std::vector<int> &island,
                           int ncount)
try : lp_vec(lp_x)
{
    integral = true;

    for (int i = 0; i < lp_vec.size(); ++i)
        if (lp_vec[i] >= Eps::Zero) {
            support_indices.push_back(i);
            support_ecap.push_back(lp_vec[i]);
            support_elist.push_back(edges[i].end[0]);
            support_elist.push_back(edges[i].end[1]);

            if (lp_vec[i] <= 1 - Eps::Zero)
                integral = false;
        }

    supp_graph = Graph::AdjList(ncount, edges, lp_x, support_indices);
    connected = supp_graph.connected(island, 0);
} catch (const exception &e) {
    cerr << e.what() << "\n";
    throw runtime_error("SupportGroup constructor failed.\n");
}

SupportGroup::SupportGroup(SupportGroup &&SG) noexcept
    : lp_vec(std::move(SG.lp_vec)),
      support_indices(std::move(SG.support_indices)),
      support_elist(std::move(SG.support_elist)),
      support_ecap(std::move(SG.support_ecap)),
      supp_graph(std::move(SG.supp_graph)),
      connected(SG.connected),
      integral(SG.integral)
{}

SupportGroup &SupportGroup::operator=(SupportGroup &&SG) noexcept
{
    lp_vec = std::move(SG.lp_vec);
    support_indices = std::move(SG.support_indices);
    support_elist = std::move(SG.support_elist);
    support_ecap = std::move(SG.support_ecap);
    supp_graph = std::move(SG.supp_graph);
    connected= SG.connected;
    integral = SG.integral;

    return *this;
}

bool SupportGroup::in_subtour_poly()
{
    if (!connected)
        return false;

    double cutval = 2.0;
    double rhs = 2.0 - Eps::MinCut;

    if (CCcut_mincut(supp_graph.node_count, support_ecap.size(),
                     &support_elist[0], &support_ecap[0], &cutval, NULL, NULL))
        throw runtime_error("CCcut_mincut failed running in_subtour_poly");

    return cutval > rhs;
}

}
}
