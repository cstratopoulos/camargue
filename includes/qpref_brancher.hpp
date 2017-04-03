/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */ /**
 * @file
 * @brief Template classes/methods for branching by priority queue.
 * @see abc_nodesel.hpp for compile-time instantiations of best tour and best
 * bound branching.
 */ /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef CMR_QPREF_BRANCHER_H
#define CMR_QPREF_BRANCHER_H

#include "datagroups.hpp"
#include "active_tour.hpp"
#include "branch_node.hpp"
#include "base_brancher.hpp"
#include "err_util.hpp"

#include <vector>
#include <queue>
#include <stdexcept>

namespace CMR {
namespace ABC {

/// Class template for branching with priority queue via some preference rule.
/// @tparam q_pref the ranking function used to instantiate branching
/// relative to a specific node selection criterion.
template <BranchNode::Pref q_pref>
class QprefBrancher : public BaseBrancher {
public:
    QprefBrancher(const Data::Instance &inst, const LP::ActiveTour &active_tour,
                  const Data::BestGroup &best_data,
                  const Graph::CoreGraph &core_graph, LP::CoreLP &core_lp);

    BranchHistory::iterator next_prob();

protected:
    void fetch_next();
    void enqueue_split(BranchNode::Split prob_array);

private:
    std::priority_queue<BranchHistory::iterator,
                        std::vector<BranchHistory::iterator>,
                        BranchNode::Pref> prob_q;
};

template<BranchNode::Pref q_pref>
QprefBrancher<q_pref>::QprefBrancher(const Data::Instance &inst,
                                     const LP::ActiveTour &active_tour,
                                     const Data::BestGroup &best_data,
                                     const Graph::CoreGraph &core_graph,
                                     LP::CoreLP &core_lp) try
    : BaseBrancher(inst, active_tour, best_data, core_graph, core_lp),
      prob_q(q_pref)
{
    prob_q.push(branch_history.begin());
} catch (const std::exception &e) {
    std::cerr << e.what() << std::endl;
    throw std::runtime_error("QprefBrancher constructor failed.");
}

template <BranchNode::Pref q_pref>
void QprefBrancher<q_pref>::enqueue_split(BranchNode::Split prob_array) try
{
    for (BranchNode &B : prob_array) {
        branch_history.emplace_front(std::move(B));
        prob_q.push(branch_history.begin());
    }
} catch (const std::exception &e) {
    std::cerr << e.what() << " putting nodes in history" << std::endl;
    throw std::runtime_error("QprefBrancher::enqueue_split failed.");
}

template <BranchNode::Pref q_pref>
BranchHistory::iterator QprefBrancher<q_pref>::next_prob()
{
    using std::cout;
    using std::endl;

    if (verbose)
        cout << "Calling QprefBrancher::next_prob...." << endl;

    if (next_itr == branch_history.end())
        fetch_next();
    else if (verbose)
        cout << "....already set." << endl;

    BranchHistory::iterator result = next_itr;

    next_itr = branch_history.end();
    return result;
}

template <BranchNode::Pref q_pref>
void QprefBrancher<q_pref>::fetch_next()
{
    using std::cout;
    using std::endl;

    if (verbose)
        cout << "Calling QprefBrancher::fetch_next..." << endl;
    if (branch_history.empty() || prob_q.empty()) {
        cout << "set to END" << endl;
        next_itr = branch_history.end();
        return;
    }

    next_itr = prob_q.top();
    prob_q.pop();
    cout << "Set next_itr to " << bnode_brief(*next_itr) << " and popped"
         << endl;
}


}
}

#endif
