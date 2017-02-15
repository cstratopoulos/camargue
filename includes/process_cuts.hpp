#ifndef CMR_PROCESS_CUTS_H
#define CMR_PROCESS_CUTS_H

#include "datagroups.hpp"
#include "cc_lpcuts.hpp"
#include "graph.hpp"
#include "tooth.hpp"
#include "cut_structs.hpp"
#include "util.hpp"

#include <memory>
#include <vector>
#include <limits>
#include <list>
#include <utility>



namespace CMR {
namespace Sep {

/** Class template for dealing with a queue of cuts of some representation. */
template<typename cut_rep>
class CutQueue {
public:
    /** Construct a CutQueue with unlimited capacity. */
    CutQueue() : q_cap(std::numeric_limits<int>::max()) {}

    /** Construct a CutQueue with capacity \p cap. */
    CutQueue(const int cap) : q_cap(cap) {}

    CutQueue(CutQueue &&CQ) noexcept = default;

    CutQueue &operator=(CutQueue &&CQ) noexcept = default;

    int q_capacity() const { return q_cap; }

    /** A reference to the most recently added cut. */
    const cut_rep &peek_front() const { return cut_q.front(); }

    /** Push a new cut to the front, popping from the back if at capacity. */
    void push_front(const cut_rep &H)
        {
            cut_q.push_front(H);
            if(cut_q.size() > q_capacity()) cut_q.pop_back();
        }

    /** Push to the back, popping from back first if at capacity. */
    void push_back(const cut_rep &H)
        {
            if(cut_q.size() >= q_capacity()) cut_q.pop_back();
            cut_q.push_back(H);
        }

    void pop_front() { cut_q.pop_front(); }  /**< Pop the front cut. */

    /** Add the cuts in Q to this list, emptying Q. */
    void splice(CutQueue<cut_rep> &Q){ cut_q.splice(cut_q.end(), Q.cut_q); }


    bool empty() const { return cut_q.empty(); } /**< Is the queue empty. */
    int size() const { return cut_q.size(); } /**< Number of cuts in queue. */

    void clear() { cut_q.clear(); } /**< Clear the queue. */

    using Itr = typename std::list<cut_rep>::iterator;
    using ConstItr = typename std::list<cut_rep>::const_iterator;

    Itr begin() { return cut_q.begin(); }
    Itr end() { return cut_q.end(); }

    ConstItr begin() const { return cut_q.begin(); }
    ConstItr end() const { return cut_q.end(); }

private:
    std::list<cut_rep> cut_q;
    int q_cap;
};

class CutTranslate {
public:
    CutTranslate(Data::GraphGroup &graph_group) :
        core_graph(graph_group.core_graph),
        delta(graph_group.delta),
        node_marks(graph_group.node_marks) {}

    void get_sparse_row(const CCtsp_lpcut_in &cc_cut,
                        const std::vector<int> &perm,
                        std::vector<int> &rmatind,
                        std::vector<double> &rmatval, char &sense,
                        double &rhs);

    void get_sparse_row(const dominoparity &dp_cut,
                        const std::vector<int> &tour_nodes,
                        std::vector<int> &rmatind, std::vector<double> &rmatval,
                        char &sense, double &rhs);

    void get_sparse_row(const std::vector<int> &handle_delta,
                        const std::vector<std::vector<int>> &tooth_edges,
                        std::vector<int> &rmatind,
                        std::vector<double> &rmatval,
                        char &sense, double &rhs);

    template<typename number_type>
    void get_activity(double &activity, const std::vector<number_type> &x,
                      const std::vector<int> &rmatind,
                      const std::vector<double> &rmatval)
        {
            activity = 0;
            for(auto i = 0; i < rmatind.size(); i++){
                int index = rmatind[i];
                activity += x[index] * rmatval[i];
            }
        }

private:
    const Graph::CoreGraph &core_graph;
    std::vector<int> &delta;
    std::vector<int> &node_marks;
};

/// Gets the indices of the teeth for an ex_blossom \p B relative to \p edges.
std::vector<int> teeth_inds(const ex_blossom &B,
                            const std::vector<int> &tour_edges,
                            const std::vector<double> &lp_vec,
                            const std::vector<Graph::Edge> &edges,
                            int ncount);

/// Like the other version, but if we already have handle_delta.
std::vector<int> teeth_inds(const ex_blossom &B,
                            const std::vector<int> &tour_edges,
                            const std::vector<double> &lp_vec,
                            const std::vector<Graph::Edge> &edges,
                            int ncount, const std::vector<int> &handle_delta);


/// Returns true if the blossom is invalid for some reason.
bool bad_blossom(const ex_blossom &B,
                 const std::vector<int> &tour_edges,
                 const std::vector<double> &lp_vec,
                 const std::vector<Graph::Edge> &edges, int ncount);

}
}



#endif
