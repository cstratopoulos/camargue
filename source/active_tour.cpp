#include "active_tour.hpp"
#include "err_util.hpp"

#include <algorithm>
#include <iostream>
#include <utility>
#include <stdexcept>
#include <string>

#include <cmath>

using std::vector;

using std::cout;
using std::cerr;
using std::endl;

using std::runtime_error;
using std::logic_error;
using std::exception;

namespace CMR {

namespace Eps = Epsilon;

namespace LP {

/**
 * Constructs an ActiveTour with basis using the Padberg-Hong algorithm for
 * getting a basis from a tour for an LP relaxation with only the degree
 * constraints.
 * @param graph the CoreGraph on which the tour was found.
 * @param best_data the data which will be used to directly initialize all
 * members except tour_base.
 */
ActiveTour::ActiveTour(const Graph::CoreGraph &graph,
                       const Data::BestGroup &best_data) try
    : tour_len(best_data.min_tour_value),
      tour_nodes(best_data.best_tour_nodes),
      perm(best_data.perm)
{
    tour_edges = vector<double>(best_data.best_tour_edges.begin(),
                                best_data.best_tour_edges.end());

    int ncount = tour_nodes.size();
    int ecount = tour_edges.size();

    vector<int> &colstat = tour_base.colstat;

    colstat = vector<int>(ecount, BStat::AtLower);
    tour_base.rowstat = vector<int>(ncount, BStat::AtLower);

    for (int i = 0; i < ncount; ++i) {
        EndPts e(tour_nodes[i], tour_nodes[(i + 1) % ncount]);
        int find_ind = graph.find_edge_ind(e.end[0], e.end[1]);

        if (find_ind == -1) {
            cerr << e << " not in core graph." << endl;
            throw logic_error ("CoreGraph does not contain all edges in tour.");
        }

        colstat[find_ind] = BStat::Basic;
    }

    if ((ncount % 2) == 0) {
        EndPts e(tour_nodes[ncount -2], tour_nodes[ncount -1]);
        int find_ind = graph.find_edge_ind(e.end[0], e.end[1]);

        colstat[find_ind] = BStat::AtUpper;

        EndPts f(tour_nodes[0], tour_nodes[ncount - 2]);
        int basis_ind = graph.find_edge_ind(f.end[0], f.end[1]);

        if (basis_ind == -1) {
            cerr << f << " not in core graph." << endl;
            throw logic_error("Graph does not contain extra basis edge.");
        }

        colstat[basis_ind] = BStat::Basic;
    }
} catch (const exception &e) {
    cerr << e.what() << endl;
    throw runtime_error("ActiveTour Padberg-Hong constructor failed.");
}

/**
 * Constructs an ActiveTour using the data that would be obtained from an
 * augmenting CoreLP::primal_pivot. Performs simple chceks to verify the tour.
 * @param tour_nodes_ this should be the "island" vector from a depth-first
 * search which verified that LP support graph was integral and connected.
 * @param lp_edges a connected, integral LP solution.
 * @param base a Basis for lp_edges.
 * @param graph the CoreGraph giving the edge set active for lp_edges, and
 * with which we verify the tour objective value and structure.
 * @remark All non-const parameters are moved from.
 */
ActiveTour::ActiveTour(vector<int> tour_nodes_, vector<double> lp_edges,
                       Basis base, double lp_objval,
                       const Graph::CoreGraph &graph) try
    : tour_nodes(std::move(tour_nodes_)), tour_edges(std::move(lp_edges)),
      tour_base(std::move(base))
{
    int ncount = graph.node_count();
    int one_count = 0;
    double manual_objval = 0.0;

    for (int i = 0; i < tour_edges.size(); ++i) {
        double &val = tour_edges[i];
        if (val > 1 - Eps::Zero) {
            ++one_count;
            manual_objval += graph.get_edge(i).len;
            val = 1.0;
        } else if (val < Eps::Zero) {
            val = 0.0;
        } else {
            cerr << "Entry " << val << " not integral." << endl;
            throw runtime_error("Updated ActiveTour with fractional vec");
        }
    }

    if (one_count != ncount) {
        cerr << one_count << " edges at one, " << ncount << " nodes." << endl;
        throw runtime_error("Wrong number of edges in ActiveTour edge vec.");
    }

    if (fabs(manual_objval - lp_objval) >= Eps::Zero) {
        cerr << "Manual objval " << manual_objval << ", lp objval "
             << lp_objval << endl;
        throw runtime_error("Objval disagreement");
    }

    tour_len = manual_objval;

    perm.resize(ncount);
    for (int i = 0; i < ncount; ++i)
        perm[tour_nodes[i]] = i;

} catch (const exception &e) {
    cerr << e.what() << endl;
    throw runtime_error("ActiveTour augmenting pivot constructor failed.");
}

/**
 * Construct an ActiveTour from scratch, using pure combinatorial data, and
 * then instate it in an LP::Relxation.
 * @param tour_nodes_ the sequence of nodes defining the tour.
 * @param relax the LP::Relaxation in which to instate the tour.
 * @param graph the CoreGraph describing the edge set of \p relax.
 * @remark \p tour_nodes_ is moved from.
 */
ActiveTour::ActiveTour(std::vector<int> tour_nodes_,
                       LP::Relaxation &relax,
                       const Graph::CoreGraph &graph) try
    : tour_len(0.0), tour_nodes(std::move(tour_nodes_))
{
    int ncount = graph.node_count();

    perm.resize(ncount);
    for (int i = 0; i < ncount; ++i)
        perm[tour_nodes[i]] = i;

    graph.tour_edge_vec(tour_nodes, tour_edges, tour_len);

    relax.copy_start(tour_edges);
    relax.factor_basis();

    double objval = relax.get_objval();
    if (fabs(objval - tour_len) >= Eps::Zero) {
        cerr << "Disagreement " << objval << " vs tour length "
             << tour_len << endl;
        throw runtime_error("ActiveTour not instated.");
    }

    tour_base = relax.basis_obj();
} catch (const exception &e) {
    cerr << e.what() << endl;
    throw runtime_error("ActiveTour new tour constructor failed.");
}

ActiveTour::ActiveTour(ActiveTour &&T) noexcept
    : tour_len(T.tour_len),
      tour_nodes(std::move(T.tour_nodes)),
      perm(std::move(T.perm)),
      tour_edges(std::move(T.tour_edges)),
      tour_base(std::move(T.tour_base))
{
    T.tour_len = 0.0;
}

ActiveTour &ActiveTour::operator=(ActiveTour &&T) noexcept
{
    tour_len = T.tour_len;
    T.tour_len = 0.0;

    tour_nodes = std::move(T.tour_nodes);
    perm = std::move(T.perm);
    tour_edges = std::move(T.tour_edges);
    tour_base = std::move(T.tour_base);

    return *this;
}

/**
 * Instate this ActiveTour in \p relax by constructing a new basis for it and
 * setting said basis (along with this tour) as the starting solution.
 */
void ActiveTour::reset_instate(Relaxation &relax)
{
    runtime_error err("Problem in ActiveTour::reset_instate");

    if (relax.num_cols() > tour_edges.size())
        tour_edges.resize(relax.num_cols(), 0.0);

    try {
        relax.copy_start(tour_edges);
        relax.factor_basis();
    } CMR_CATCH_PRINT_THROW("copying and factoring", err);

    double objval = relax.get_objval();
    if (fabs(objval - tour_len) >= Eps::Zero) {
        cerr << "Tour len " << tour_len << " vs lp objval "
             << objval << endl;
        throw err;
    }

    try { tour_base = relax.basis_obj(); }
    CMR_CATCH_PRINT_THROW("getting the basis", err);
}

/**
 * This function will try to instate this tour in \p relax, first by copying
 * the tour edges and basic statuses, and then by generating a new basis
 * if that fails.
 * @post the tour is instated.
 */
void ActiveTour::instate(Relaxation &relax)
{
    runtime_error err("Problem in ActiveTour::instate");

    try {
        relax.copy_start(tour_edges, tour_base.colstat, tour_base.rowstat);
        relax.factor_basis();
    } CMR_CATCH_PRINT_THROW("copying full start and factoring", err)

    double objval = relax.get_objval();
    if (fabs(objval - tour_len) < Eps::Zero)
        return;

    try { reset_instate(relax); }
    CMR_CATCH_PRINT_THROW("doing reset after full copy_start failed", err);
}

/**
 * This method will update the BestGroup \p best_data with the stored
 * information in this ActiveTour. This method will verify that the tour
 * is indeed an improvement on the one currently in \p best_data, throwing
 * an exception if it is worse or if the edge vector is generated incorrectly.
 */
void ActiveTour::best_update(Data::BestGroup &best_data) const
{
    if (tour_len > best_data.min_tour_value)
        throw runtime_error("Tried AciveTour::best_update with worse tour!!");

    runtime_error err("Problem in ActiveTour::best_update");

    best_data.min_tour_value = tour_len;
    try {
        best_data.best_tour_nodes = tour_nodes;
        best_data.perm = perm;

        best_data.best_tour_edges = vector<int>(tour_edges.begin(),
                                                tour_edges.end());
    } CMR_CATCH_PRINT_THROW("copying vectors", err);

    int onecount = std::count(best_data.best_tour_edges.begin(),
                              best_data.best_tour_edges.end(), 1);
    int ncount = tour_nodes.size();

    if (ncount != onecount) {
        std::string mismatch(std::to_string(ncount) + " nodes but "
                             + std::to_string(onecount) +
                             " edges equal one in ActiveTour::best_update");
        throw runtime_error(mismatch);
    }
}

}
}
