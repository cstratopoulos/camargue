#include "tooth.hpp"
#include "config.hpp"
#include "err_util.hpp"

#if CMR_HAVE_TIMSORT
#include <timsort.hpp>
#endif

extern "C" {
#include <concorde/INCLUDE/cut.h>
}

#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <tuple>

using std::array;
using std::vector;

using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::to_string;

using std::pair;
using std::unique_ptr;

using std::runtime_error;
using std::exception;

namespace CMR {
namespace Sep {

using ToothList = CandidateTeeth::ToothList;
using IteratorMat = CandidateTeeth::IteratorMat;

static inline bool ptr_cmp(const SimpleTooth::Ptr &S, const SimpleTooth::Ptr &T)
{ return S->body_size() < T->body_size(); }

static inline bool elim_less_tie(const SimpleTooth::Ptr &S,
				 const SimpleTooth::Ptr &T)
{
    return std::make_tuple(S->slack, S->body_size()) <
        std::make_tuple(T->slack, T->body_size());
}

static void tlist_sort(ToothList &T)
{
#if CMR_HAVE_TIMSORT
    gfx::timsort(T.begin(), T.end(), ptr_cmp);
#else
    std::sort(T.begin(), T.end(), ptr_cmp);
#endif
}

vector<vector<int>> CandidateTeeth::adj_zones;
vector<IteratorMat> CandidateTeeth::seen_ranges;

CandidateTeeth::CandidateTeeth(const LP::ActiveTour &active_tour_,
			       Data::SupportGroup &_supp_dat) :
    light_teeth(std::vector<ToothList>(_supp_dat.supp_graph.node_count)),
    list_sizes(_supp_dat.supp_graph.node_count, {{0, 0, 0}}),
    endmark(_supp_dat.supp_graph.node_count, CC_LINSUB_BOTH_END),
    active_tour(active_tour_),
    supp_dat(_supp_dat),
    t_all("Candidate Teeth"),
    t_zones("Adj zones", &t_all),
    t_find("Initial find", &t_all),
    t_sort("Sort by root", &t_all)
{
    t_all.start();
    t_zones.start();
    Graph::AdjList &G_s = supp_dat.supp_graph;
    int ncount = G_s.node_count;
    const vector<int> &perm = active_tour.tour_perm();
    const vector<int> &tour = active_tour.nodes();

    seen_ranges.resize(ncount);
    adj_zones.resize(ncount, vector<int>(ncount, 0));

    std::fill(adj_zones.begin(), adj_zones.end(),
              vector<int>(ncount, 0));

#ifdef CMR_USE_OMP
#pragma omp parallel for
#endif
    for (int root_ind = 0; root_ind < ncount; ++root_ind) {
        int actual_vx = tour[root_ind];
        const Graph::Node &x = G_s.nodelist[actual_vx];

        // this needs to be done for iterator validity
        light_teeth[root_ind].reserve(2 * (x.degree() - 1));

        for (const Graph::AdjObj &a : x.neighbors) {
            int end_ind = perm[a.other_end];
            adj_zones[root_ind][end_ind] = 1;
        }

        int label = 0;
        for (int &i : adj_zones[root_ind]) {
            if (i == 1) {
                ++label;
                i = -1 * label;
            } else
                i = label;
        }

        seen_ranges[root_ind] = IteratorMat(label + 1,
                                            light_teeth[root_ind].rend());
    }

    t_zones.stop();
}

void CandidateTeeth::get_light_teeth()
{
    t_find.start();
    runtime_error err("Problem in CandidateTeeth::get_light_teeth.");
    unique_ptr<LinsubCBData> cb_data;
    Graph::AdjList &supp_graph = supp_dat.supp_graph;
    vector<bool> node_marks(supp_graph.node_count, false);

    const vector<int> &tour_nodes = active_tour.nodes();
    const vector<int> &perm = active_tour.tour_perm();


    try {
        cb_data =
        util::make_unique<LinsubCBData>(light_teeth,
                                        adj_zones, seen_ranges,
                                        list_sizes,
                                        node_marks,
                                        tour_nodes, perm,
                                        supp_graph);
    } CMR_CATCH_PRINT_THROW("allocating LinsubCBData.", err);

    if (CCcut_linsub_allcuts(supp_graph.node_count, supp_graph.edge_count,
                             const_cast<int *>(&tour_nodes[0]), &endmark[0],
                             &supp_dat.support_elist[0],
                             &supp_dat.support_ecap[0],
                             3.0 - Epsilon::Cut,
                             cb_data.get(), teeth_cb))
        throw err;

    t_find.stop();
}

void CandidateTeeth::sort_by_root()
{
    t_sort.start();

#ifdef CMR_USE_OMP
#pragma omp parallel for
#endif
    for (int root = 0; root < light_teeth.size(); ++root) {
        ToothList &teeth = light_teeth[root];
        if (teeth.empty())
            continue;

        array<int, 3> &sizes = list_sizes[root];

        bool have_left = (sizes[0] >= 0);
        bool have_right = (sizes[1] >= 0);
        bool have_dist = (sizes[2] >= 0);
        bool do_sort = (have_left + have_right + have_dist >= 2);

        if (do_sort)
            tlist_sort(teeth);
    }

    t_sort.stop();
    t_all.stop();
}

IntPair CandidateTeeth::get_range(const int root, const ToothBody &s,
                                  const vector<vector<int>> &zones)
{
    IntPair range = IntPair(-1, -1);
    int start = zones[root][s.start], end = zones[root][s.end];

    if (start == end) {//if same endpoints
        if (start < 0) {//singleton
            range = IntPair(-start, -start);
        }//else empty
        return range;
    }

    //now different endpoints
    if (start < 0) {
        range.first = -start; //first definitely start
        range.second = fabs(end);
        return range;
    }

    //now diff endpoints, start >= 0
    if (end < 0) {
        range.second = -end; //second definitely end
        if (start == -end)
            range.first = start; //start is itself
        else
            range.first = start + 1;// fmin(start + 1, -end) start is at most end
        return range;
    }

    //now diff endpoints, both >= 0
    range = IntPair(start + 1, end);// fmax(start + 1, end - 1));
    return range;
}

bool CandidateTeeth::root_equivalent(const int root, const ToothBody &s1,
				     const ToothBody &s2,
				     const vector<vector<int>> &zones)
{
    return get_range(root, s1, zones) == get_range(root, s2, zones);
}

bool CandidateTeeth::root_equivalent(const int root, const ToothBody &s1,
				     const ToothBody &s2) const
{
    return root_equivalent(root, s1, s2, adj_zones);
}

int CandidateTeeth::teeth_cb(double cut_val, int cut_start, int cut_end,
			     void *u_data)
{

    int rval = 0;

    double slack = (cut_val - (2.0 - Epsilon::Cut)) / 2.0;

    LinsubCBData *arg = (LinsubCBData *) u_data;
    if ((cut_end - cut_start + 1) > (arg->G_s.node_count - 2))
        return 0;

    vector<bool> &marks = arg->node_marks;
    const vector<int> &tour = arg->tour_nodes;
    const vector<int> &perm = arg->perm;

    const Graph::AdjList &G = arg->G_s;
    int ncount = G.node_count;

    std::unordered_map<int, double> &rb_sums = arg->rb_sums;
    double rb_lower = cut_val - (1.5 - Epsilon::Cut);

    ToothBody &old_seg = arg->old_seg;
    vector<vector<int>> &zones = arg->adj_zones;
    vector<IteratorMat> &ranges = arg->ranges;

    //distant add
    if (cut_start == old_seg.start) {
        for (int i = old_seg.end + 1; i <= cut_end; ++i) {
            marks[i] = true;
            rb_sums.erase(i);
        }

        for (int i = old_seg.end + 1; i <= cut_end; ++i) {
            const Graph::Node &vx = G.nodelist[tour[i]];
            for (const Graph::AdjObj &a : vx.neighbors) {
                int root_perm = perm[a.other_end];

                if (marks[root_perm] == false)
                    rb_sums[root_perm] += a.val;
            }
        }
    } else {//on a new degree node
        //clean up after the old one
        std::fill(marks.begin(), marks.end(), false);
        rb_sums.clear();

        //set up for the new one
        marks[cut_start] = true;

        const Graph::Node &vx = G.nodelist[tour[cut_start]];
        for (const Graph::AdjObj &a : vx.neighbors) {
            int root_perm = perm[a.other_end];

            if (marks[root_perm] == false)
                rb_sums[root_perm] += a.val;
        }
    }

    for (auto &kv : rb_sums) {
        int root = kv.first;
        double rb_sum = kv.second;
        if ((cut_start == cut_end) && (root > cut_end) && (root != (ncount - 1)))
            continue;
        if (((cut_end - cut_start + 1) == (ncount - 2)) && (root > cut_end))
            continue;

        if (rb_sum > rb_lower) {
            ToothList &teeth = arg->light_teeth[root];
            double abs_slack = fabs(cut_val - rb_sum - 1);
            double new_slack = (abs_slack < Epsilon::Zero) ? 0 : abs_slack;

            try {
                add_tooth(teeth, zones, ranges, arg->list_sizes[root],
                          root, cut_start, cut_end, new_slack);
            } catch (const exception &e) {
                cerr << e.what() << " pushing back tooth in teeth_cb.\n";
                rval = 1;
                break;
            }
        }
    }

    old_seg = ToothBody(cut_start, cut_end, slack); //right adjacent update
    return rval;
}

inline void CandidateTeeth::add_tooth(ToothList &teeth,
                                      const vector<vector<int>> &zones,
                                      vector<IteratorMat> &ranges,
                                      array<int, 3> &sizes,
                                      const int root, const int body_start,
                                      const int body_end, const double slack)
{
    bool elim = false;
    ToothBody body(body_start, body_end);
    IntPair range = get_range(root, body, zones);

    if (!teeth.empty()) {
        ToothList::reverse_iterator &rit = ranges[root](range.first,
                                                        range.second);

        if (rit != teeth.rend()) {
            elim = true;
            double old_slack = (*rit)->slack;
            if (slack < old_slack ||
                (slack <= old_slack && body.size() < (*rit)->body_size())) {
                --sizes[static_cast<int>((*rit)->type())];
                *rit = util::make_unique<SimpleTooth>(root, body, slack);
                ++sizes[static_cast<int>((*rit)->type())];
            }
        }
    }

    if (!elim) {
        teeth.emplace_back(util::make_unique<SimpleTooth>(root, body, slack));
        ranges[root](range.first, range.second) = teeth.rbegin();
        ++sizes[static_cast<int>(teeth.back()->type())];
    }
}

string CandidateTeeth::print_label(const SimpleTooth &T)
{
    string
    body = (T.body_start == T.body_end) ?
    to_string(T.body_start) :
    ("{"
     + to_string(T.body_start) + "..." + to_string(T.body_end)
     + "}");

    return "(" + to_string(T.root) + ", " + body + ")";
}

void CandidateTeeth::print_tooth(const SimpleTooth &T, bool full)
{
    print_tooth(T, full, active_tour.nodes());
}

void CandidateTeeth::print_tooth(const SimpleTooth &T, bool full,
				 const vector<int> &bt)
{
    cout << "(" << bt[T.root] << ", {" << bt[T.body_start];

    if (T.body_start == T.body_end) {
        cout << "}) -- slack " << T.slack << "\n";
        return;
    }

    if (!full) {
        cout << ", ..., " << bt[T.body_end] << "}) -- slack " << T.slack << "\n";
        return;
    }

    bool comma_sep = T.body_end - T.body_start <= 20;
    int i = T.body_start;

    cout << ", ";
    while (i++ != T.body_end) {
        cout << bt[i];
        if (i == T.body_end)
            cout << "}";
        else
            cout << (comma_sep ? ", " : "\n\t");
    }
    cout << ") -- slack " << T.slack << "\n";
}

void CandidateTeeth::profile()
{
    t_zones.report(true);
    t_find.report(true);
    t_sort.report(true);
    t_all.report(true);
}

}
}
