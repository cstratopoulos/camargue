#ifndef CMR_TOOTH_HPP
#define CMR_TOOTH_HPP

#include "Graph.hpp"
#include "datagroups.hpp"
#include "util.hpp"
#include "timer.hpp"

#include <memory>
#include <vector>
#include <string>
#include <utility>
#include <unordered_map>
#include <bitset>

namespace CMR {

enum class ListStat {None, Merge, Full};

struct ToothBody : Segment {
    ToothBody() = default;
    ToothBody(int _start, int _end, double _slack) :
        Segment(_start, _end), slack(_slack) {}
    ToothBody(int _start, int _end) : Segment(_start, _end), slack(1.0) {}

    double slack;
};

struct tooth_seg {
  tooth_seg(int _start, int _end, double _slack) :
    start(_start), end(_end), slack(_slack) {}
  tooth_seg(int _start, int _end) : start(_start), end(_end), slack(1.0) {}

  int body_size() const { return end - start + 1; }
  bool contains(int vx) const { return start <= vx && vx <= end; }
  bool subset_of(tooth_seg &S) const {
    return S.start <= start && end <= S.end;
  }

  int start, end;
  double slack;
};

struct SimpleTooth {
  SimpleTooth(int _root, int _body_start, int _body_end, double _slack) :
    root(_root), body_start(_body_start), body_end(_body_end),
    slack(_slack) {}

  SimpleTooth(int _root, tooth_seg &seg, double _slack) :
    root(_root), body_start(seg.start), body_end(seg.end), slack(_slack) {}

  typedef std::unique_ptr<SimpleTooth> Ptr;

  int root, body_start, body_end;
  int cutgraph_index;
  double slack;

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
  CandidateTeeth(CMR::Data::GraphGroup &_graph_dat,
		 CMR::Data::BestGroup &_best_dat,
		 CMR::Data::SupportGroup &_supp_dat);

  int get_light_teeth();

  static void get_range(const int root, const tooth_seg &s, IntPair &range,
			const std::vector<std::vector<int>> &zones);

  static bool root_equivalent(const int root, const tooth_seg &s1,
			      const tooth_seg &s2,
			      const std::vector<std::vector<int>> &zones);

  bool root_equivalent(const int root, const tooth_seg &s1,
		       const tooth_seg &s2) const;

  int merge_and_sort(const int root);
  int merge_and_sort();
  
  void unmerged_weak_elim();
  void complement_elim();
  
  static void print_tooth(const SimpleTooth &T, bool full,
			  const std::vector<int> &tour_nodes);
  void print_tooth(const SimpleTooth &T, bool full);
  std::string print_label(const SimpleTooth &T);

  void profile();

  std::vector<std::vector<SimpleTooth::Ptr>>
  light_teeth, left_teeth, right_teeth, dist_teeth;

  static std::vector<std::vector<int>> adj_zones;
  std::vector<ListStat> stats;

private:
  friend class DPCutGraph;
  friend class DPwitness;
  
  std::vector<int> endmark;
  
  static void add_tooth(std::vector<SimpleTooth::Ptr> &teeth,
			const std::vector<std::vector<int>> &zones,
			const int root, const int body_start,
			const int body_end, const double slack);
  
  static int teeth_cb(double cut_val, int cut_start, int cut_end,
		      void *u_data);

  struct LinsubCBData {
    LinsubCBData(std::vector<std::vector<SimpleTooth::Ptr>> &_r_teeth,
		 std::vector<std::vector<SimpleTooth::Ptr>> &_l_teeth,
		 std::vector<std::vector<SimpleTooth::Ptr>> &_d_teeth,
		 std::vector<std::vector<int>> &_adj_zones,
		 std::vector<int> &_node_marks,
		 std::vector<int> &_tour_nodes,
		 std::vector<int> &_perm,
		 CMR::SupportGraph &_G_s) :
      r_teeth(_r_teeth), l_teeth(_l_teeth), d_teeth(_d_teeth),
      adj_zones(_adj_zones),
      node_marks(_node_marks),
      tour_nodes(_tour_nodes), perm(_perm),
      G_s(_G_s),
      old_seg(_G_s.node_count - 1, _G_s.node_count - 1, 0.0),
      prev_slacks(std::vector<
		  std::pair<int, double>
		  >(_G_s.node_count,
		    std::pair<int, double>(_G_s.node_count + 1, 1.0)))
    {}

    std::vector<std::vector<SimpleTooth::Ptr>> &r_teeth, &l_teeth, &d_teeth;
    std::vector<std::vector<int>> &adj_zones;

    std::vector<int> &node_marks;
    std::vector<int> &tour_nodes, &perm;

    CMR::SupportGraph &G_s;

    tooth_seg old_seg;
    
    std::unordered_map<int, double> rb_sums;
    std::vector<std::pair<int, double>> prev_slacks;
  };

  CMR::Data::GraphGroup &graph_dat;
  CMR::Data::BestGroup &best_dat;
  CMR::Data::SupportGroup &supp_dat;

  CMR::Timer t_all;
  CMR::Timer t_zones;
  CMR::Timer t_find;
  CMR::Timer t_elim;
  CMR::Timer t_comp_elim;
  CMR::Timer t_sort;
};

}

#endif
