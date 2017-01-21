#include "brancher.hpp"
#include "err_util.hpp"
#include "util.hpp"

#include <stdexcept>

using std::vector;

using std::cout;
using std::cerr;
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

ostream &operator<<(ostream &os, Problem &prob)
{
    os << prob.status << " " << prob.type << " branch on edge "
       << prob.edge_ind;
    return os;
}

Brancher::Brancher(LP::Relaxation &lp_rel,
                   const vector<Edge> &edges, const LP::TourBasis &tbase,
                   const double &tourlen) try
    : lp_relax(lp_rel), core_edges(edges), tour_base(tbase), tour_len(tourlen)
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

        //enact_top();

        if (top.status == Pstat::Unseen)
            return top;

        subprobs.pop();
        
        if (top.status == Pstat::Seen)
            split_prob(branch_edge());
    }

    return solved_prob;
}

int Brancher::branch_edge()
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

    return winner.index;
}

void Brancher::split_prob(int edge)
{
    subprobs.emplace(Ptype::Contra, edge);
    subprobs.emplace(Ptype::Affirm, edge);
}

void Brancher::enact_top()
{
    // if (subprobs.empty() || subprobs.top().type == Ptype::Root)
    //     return;

    // const Problem &top = subprobs.top();
    // const int ind = top.edge_ind;

    // if (top.type == Ptype::Affirm) {
    //     double tour_entry = tour_base.best_tour_edges[ind];
        
    //     if (top.status == Pstat::Unseen)
    //         lp_relax.tighten_bound(ind, (tour_entry == 1) ? 'L' : 'U',
    //                                tour_entry);
    //     else 
    //         lp_relax.tighten_bound(ind, (tour_entry == 1) ? 'U' : 'L',
    //                                tour_entry);
    // } else {
    //     cout << "Not yet set up for contra branching.\n";
    // }
}

}
}
