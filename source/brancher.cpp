#include "brancher.hpp"
#include "err_util.hpp"
#include "util.hpp"

#include <stdexcept>

using std::vector;

using std::cout;
using std::cerr;
using std::endl;
using std::ostream;

using std::runtime_error;
using std::logic_error;
using std::exception;

namespace CMR {
namespace ABC {

using Edge = Graph::Edge;

using Ptype = Problem::Type;
using Pstat = Problem::Status;

Problem Brancher::solved_prob(Ptype::Root, Pstat::Pruned, -1);

Brancher::Brancher(LP::Relaxation &lp_rel,
                   const vector<Edge> &edges, const LP::TourBasis &tbase,
                   const double &tourlen, const ContraStrat strat) try
    : lp_relax(lp_rel), core_edges(edges), tour_base(tbase), tour_len(tourlen),
      contra_strategy(strat),
      contra_enforce(contra_fix_enforce), contra_undo(contra_fix_undo)
{
    subprobs.emplace(Ptype::Root, Pstat::Seen, -1);
} catch (const exception &e) {
    cerr << e.what() << "\n";
    throw runtime_error("Brancher constructor failed.");
}

Problem &Brancher::next_prob()
{
    runtime_error err("Problem in Brancher::next_prob.");

    while (!subprobs.empty()) {
        Problem &top = subprobs.top();

        enact_top();

        if (top.status == Pstat::Unseen)
            return top;

        subprobs.pop();
        
        if (top.status == Pstat::Seen)
            split_prob(branch_edge_index());
    }

    return solved_prob;
}

int Brancher::branch_edge_index()
{
    vector<double> x = lp_relax.lp_vec();
    vector<int> colstat = lp_relax.col_stat();
    vector<int> lw_inds;

    for (int i = 0; i < colstat.size(); ++i)
        if (colstat[i] == 1 && !util::var_integral(x[i]))
            lw_inds.push_back(i);

    vector<int> sb1inds = length_weighted_cands(core_edges, lw_inds, x,
                                                SB1Cands);
    vector<double> downobj;
    vector<double> upobj;

    lp_relax.primal_strong_branch(tour_base.best_tour_edges,
                                  tour_base.colstat, tour_base.rowstat,
                                  sb1inds, downobj, upobj, SB1Lim, tour_len);

    vector<ScoreTuple> sb2cands = ranked_cands(sb1inds, downobj, upobj,
                                               StrongMult, tour_len, SB2Cands);
    vector<int> sb2inds;

    for (ScoreTuple &t : sb2cands)
        sb2inds.push_back(t.index);

    lp_relax.primal_strong_branch(tour_base.best_tour_edges,
                                  tour_base.colstat, tour_base.rowstat,
                                  sb2inds, downobj, upobj, SB2Lim, tour_len);

    ScoreTuple winner = ranked_cands(sb2inds, downobj, upobj, StrongMult,
                                     tour_len, 1)[0];
    cout << "\tWinner estimates: " << winner.down_est << ", "
         << winner.up_est << "\n";
    return winner.index;
}

void Brancher::split_prob(int edge)
{
    subprobs.emplace(Ptype::Contra, edge);
    subprobs.emplace(Ptype::Affirm, edge);
}

void Brancher::enact_top()
{
    if (subprobs.empty() || subprobs.top().type == Ptype::Root)
        return;

    const Problem &top = subprobs.top();
    const int ind = top.edge_ind;
    const double tour_entry = tour_base.best_tour_edges[ind];

    if (top.type == Ptype::Affirm) {        
        if (top.status == Pstat::Unseen) { // Clamp the edge to affirm tour.
            cout << "Clamping " << top << endl;
            if (tour_entry == 1)
                lp_relax.tighten_bound(ind, 'L', 1);
            else if (tour_entry == 0)
                lp_relax.tighten_bound(ind, 'U', 0);
        } else { //Undo the clamp.
            cout << "Unclamping " << top << endl;
            if (tour_entry == 1)
                lp_relax.tighten_bound(ind, 'L', 0);
            else if (tour_entry == 0)
                lp_relax.tighten_bound(ind, 'U', 1);
        }
    } else {
        if (top.status == Pstat::Unseen) {
            cout << "Enforcing " << top << endl;
            contra_enforce(lp_relax, ind, tour_entry);
        } else {
            cout << "Unenforcing " << top << endl;
            contra_undo(lp_relax, ind, tour_entry);
        }
    }
}

/**
 * @param[in] rel the Relaxation to modify.
 * @param[in] branch_ind the index of the edge to branch on.
 * @param[in] tour_entry the entry of the best tour corresponding to 
 * \p branch_ind.
 * This function modifies the Relaxation to explore a branch where the solution
 * is fixed to disagree with the current best tour. This is done by changing
 * bounds on \p branch_ind so that the upper bound is zero if \p tour_entry
 * is one, and the lower bound is one if \p tour_entry is zero. In either case,
 * effectively adds the constraint that `x[branch_ind] = 1 - tour_entry` to 
 * \p rel.
 */
void contra_fix_enforce(LP::Relaxation &rel, const int branch_ind,
                             const double tour_entry)
{
    if (tour_entry == 1)
        rel.tighten_bound(branch_ind, 'U', 0);
    else if (tour_entry == 0)
        rel.tighten_bound(branch_ind, 'L', 1);
    
    return;
}

/**
 * @param[in] rel the Relaxation to modify.
 * @param[in] branch_ind the index of the edge to undo a branch on.
 * @param[in] tour_entry the entry of the best tour corresponding to 
 * \p branch_ind.
 * This function undoes a branch which was performed by contra_fix_enforce. It
 * assumes that upper or lower bounds have been modified precisely in the way
 * indicated by that function, hence restores them to normal values. If 
 * \p tour_entry is one, then branch_ind will have its upper bound relaxed
 * back to one. If it is zero, it will have its lower bound relaxed back to
 * zero.
 */
void contra_fix_undo(LP::Relaxation &rel, const int branch_ind,
                             const double tour_entry)
{
    if (tour_entry == 1)
        rel.tighten_bound(branch_ind, 'U', 1);
    else if (tour_entry == 0)
        rel.tighten_bound(branch_ind, 'L', 0);

    return;
}

ostream &operator<<(ostream &os, Ptype type)
{
    switch (type) {
    case Ptype::Root:
        os << "Root";
        break;
    case Ptype::Affirm:
        os << "Affirm";
        break;
    case Ptype::Contra:
        os << "Contra";
        break;
    }

    return os;
}

ostream &operator<<(ostream &os, Pstat stat)
{
    switch (stat) {
    case Pstat::Pruned:
        os << "Pruned";
        break;
    case Pstat::Unseen:
        os << "Unseen";
        break;
    case Pstat::Seen:
        os << "Seen";
        break;
    }
    return os;
}

ostream &operator<<(ostream &os, const Problem &prob)
{
    os << prob.status << " " << prob.type << " branch on edge "
       << prob.edge_ind;
    return os;
}

}
}
