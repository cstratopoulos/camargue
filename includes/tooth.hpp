#ifndef CMR_TOOTH_HPP
#define CMR_TOOTH_HPP

#include "graph.hpp"
#include "datagroups.hpp"
#include "cut_structs.hpp"
#include "util.hpp"
#include "timer.hpp"

#include <array>
#include <memory>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace CMR {
namespace Sep {

enum class ListStat {None, Merge, Full};

class CandidateTeeth {
public:
    CandidateTeeth(Data::BestGroup &_best_dat,
                   Data::SupportGroup &_supp_dat);

    void get_light_teeth();

    static IntPair get_range(const int root, const ToothBody &s,
                             const std::vector<std::vector<int>> &zones);

    static bool root_equivalent(const int root, const ToothBody &s1,
                                const ToothBody &s2,
                                const std::vector<std::vector<int>> &zones);

    bool root_equivalent(const int root, const ToothBody &s1,
                         const ToothBody &s2) const;

    void sort_by_root();

    static void print_tooth(const SimpleTooth &T, bool full,
                            const std::vector<int> &tour_nodes);
    void print_tooth(const SimpleTooth &T, bool full);
    std::string print_label(const SimpleTooth &T);

    void profile();

    using ToothList = std::vector<SimpleTooth::Ptr>;
    using IteratorMat = util::SquareUT<ToothList::reverse_iterator>;

    std::vector<ToothList> light_teeth;

    static std::vector<std::vector<int>> adj_zones;
    static std::vector<IteratorMat> seen_ranges;
    std::vector<std::array<int, 3>> list_sizes;
    std::vector<ListStat> stats;

private:
    friend class DPwitness;

    std::vector<int> endmark;

    static void add_tooth(ToothList &teeth,
                          const std::vector<std::vector<int>> &zones,
                          std::vector<IteratorMat> &ranges,
                          std::array<int, 3> &sizes,
                          const int root, const int body_start,
                          const int body_end, const double slack);

    static int teeth_cb(double cut_val, int cut_start, int cut_end,
                        void *u_data);

    struct LinsubCBData {
        LinsubCBData(std::vector<ToothList> &_light_teeth,
                     std::vector<std::vector<int>> &_adj_zones,
                     std::vector<IteratorMat>
                     &_ranges,
                     std::vector<std::array<int, 3>> &_list_sizes,
                     std::vector<bool> &_node_marks,
                     std::vector<int> &_tour_nodes,
                     std::vector<int> &_perm,
                     const Graph::AdjList &_G_s) :
            light_teeth(_light_teeth),
            adj_zones(_adj_zones), ranges(_ranges), list_sizes(_list_sizes),
            node_marks(_node_marks),
            tour_nodes(_tour_nodes), perm(_perm),
            G_s(_G_s),
            old_seg(_G_s.node_count - 1, _G_s.node_count - 1, 0.0)
            {}

        std::vector<std::vector<SimpleTooth::Ptr>> &light_teeth;
        std::vector<std::vector<int>> &adj_zones;
        std::vector<IteratorMat> &ranges;
        std::vector<std::array<int, 3>> &list_sizes;

        std::vector<bool> &node_marks;
        std::vector<int> &tour_nodes;
        std::vector<int> &perm;

        const Graph::AdjList &G_s;

        ToothBody old_seg;

        std::unordered_map<int, double> rb_sums;
    };

    Data::BestGroup &best_dat;
    Data::SupportGroup &supp_dat;

    CMR::Timer t_all;
    CMR::Timer t_zones;
    CMR::Timer t_find;
    CMR::Timer t_sort;
};

}
}

#endif
