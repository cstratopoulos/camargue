#ifndef CMR_TOOTH_HPP
#define CMR_TOOTH_HPP

#include "graph.hpp"
#include "datagroups.hpp"
#include "util.hpp"
#include "timer.hpp"

#include <array>
#include <memory>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace CMR {

enum class ListStat {None, Merge, Full};

struct ToothBody : Segment {
    ToothBody() = default;
    ToothBody(int _start, int _end, double _slack) :
        Segment(_start, _end), slack(_slack) {}
    ToothBody(int _start, int _end) : Segment(_start, _end), slack(1.0) {}

    double slack;
};

struct SimpleTooth {
    SimpleTooth(int _root, int _body_start, int _body_end, double _slack) :
        root(_root), body_start(_body_start), body_end(_body_end),
        slack(_slack) {}

    SimpleTooth(int _root, ToothBody &seg, double _slack) :
        root(_root), body_start(seg.start), body_end(seg.end), slack(_slack) {}

    typedef std::unique_ptr<SimpleTooth> Ptr;

    int root, body_start, body_end;
    int cutgraph_index;
    double slack;

    enum Type {
        LeftAdj = 0,
        RightAdj = 1,
        Dist = 2
    };

    Type type() const
        {
            if (body_start == root + 1)
                return LeftAdj;
            if (body_end + 1 == root)
                return RightAdj;
            return Dist;
        }

    int body_size() const { return body_end - body_start + 1; }
    bool body_contains(int i) const { return body_start <= i && i <= body_end; }
    bool is_subset_of(const SimpleTooth &T) const {
        return root == T.root &&
        T.body_start <= body_start &&
        body_end <= T.body_end;
    }
    
};

class CandidateTeeth {
public:
    CandidateTeeth(Data::GraphGroup &_graph_dat,
                   Data::BestGroup &_best_dat,
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

    std::vector<ToothList> light_teeth;

    static std::vector<std::vector<int>> adj_zones;
    static std::vector<util::SquareUT<ToothList::reverse_iterator>> seen_ranges;
    std::vector<ListStat> stats;

private:
    friend class DPwitness;
  
    std::vector<int> endmark;
    std::vector<std::array<int, 3>> list_sizes;
  
    static void add_tooth(ToothList &teeth,
                          const std::vector<std::vector<int>> &zones,
                          std::vector<
                          util::SquareUT<ToothList::reverse_iterator>
                          > &ranges,
                          std::array<int, 3> &sizes,
                          const int root, const int body_start,
                          const int body_end, const double slack);
  
  static int teeth_cb(double cut_val, int cut_start, int cut_end,
		      void *u_data);

  struct LinsubCBData {
      LinsubCBData(std::vector<ToothList> &_light_teeth,
                   std::vector<std::vector<int>> &_adj_zones,
                   std::vector<util::SquareUT<ToothList::reverse_iterator>>
                   &_ranges,
                   std::vector<std::array<int, 3>> &_list_sizes,
                   std::vector<int> &_node_marks,
                   std::vector<int> &_tour_nodes,
                   std::vector<int> &_perm,
                   CMR::SupportGraph &_G_s) :
          light_teeth(_light_teeth),
          adj_zones(_adj_zones), ranges(_ranges), list_sizes(_list_sizes),
          node_marks(_node_marks),
          tour_nodes(_tour_nodes), perm(_perm),
          G_s(_G_s),
          old_seg(_G_s.node_count - 1, _G_s.node_count - 1, 0.0)
          {}

      std::vector<std::vector<SimpleTooth::Ptr>> &light_teeth;
      std::vector<std::vector<int>> &adj_zones;
      std::vector<util::SquareUT<ToothList::reverse_iterator>> &ranges;
      std::vector<std::array<int, 3>> &list_sizes;

      std::vector<int> &node_marks;
      std::vector<int> &tour_nodes, &perm;

      CMR::SupportGraph &G_s;

      ToothBody old_seg;
    
      std::unordered_map<int, double> rb_sums;
  };

  Data::GraphGroup &graph_dat;
  Data::BestGroup &best_dat;
  Data::SupportGroup &supp_dat;

  CMR::Timer t_all;
  CMR::Timer t_zones;
  CMR::Timer t_find;
  CMR::Timer t_sort;
};

}

#endif
