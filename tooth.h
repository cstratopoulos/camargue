#ifndef PSEP_TOOTH_H
#define PSEP_TOOTH_H

#include<memory>
#include<list>
#include<algorithm>

#include "Graph.h"

class PSEP_CandTooth {
 public:
 PSEP_CandTooth(std::vector<int> & _tour_nodes, SupportGraph & _G,
		std::vector<int> & _marks) :
  best_tour_nodes(_tour_nodes), edge_marks(_marks) {
    SimpleTooth::ncount = _tour_nodes.size();
    SimpleTooth::G_s = &_G;
    SimpleTooth::best_tour_nodes = &_tour_nodes[0];

    light_teeth.resize(_tour_nodes.size());
    heavy_teeth.resize(_tour_nodes.size());
  }

  void build_collection();
  void find_root_adjacent_teeth(const int root);
  void find_root_distant_teeth(const int root);

 private:
  class SimpleTooth{
  public:
    SimpleTooth() {root = -1; body_start = -1; body_end = -1;}
  SimpleTooth(const int _root, const int _body_start, const int _body_end) :
    root(_root), body_start(_body_start), body_end(_body_end) {
      if(body_start <= body_end) // [____<---*--->__] gives sandwich
	sandwich = (body_start <= root && root <= body_end);
      else //        [--->____<--*--]  OR [--*->____<--] gives sandwich
	sandwich = (body_start <= root || root <= body_end);

      int add_one = (int) (!sandwich);
      body_size = (body_start <= body_end) ? (body_end - body_start + add_one) :
	((ncount - body_start) + body_end + add_one);
    }

    bool body_contains(const int node_perm);
    static bool C_body_subset(const SimpleTooth &T, const SimpleTooth &R);
    void complement();
    void increment_slack(const int new_vx, double *lhs_p, int *rhs_p);
    
    bool operator>(const SimpleTooth& T) const {
      return body_size > T.body_size;
    }

    static bool p_greater(const std::unique_ptr<SimpleTooth>  &T,
			  const std::unique_ptr<SimpleTooth> &R){
      return *T > *R;
    }

    void print();

  
  private:
    int root;
    int body_start;
    int body_end;
    int body_size;
    double slack;
    bool sandwich;

    friend class PSEP_CandTooth;

    //int node_index;
    static int ncount;
    static SupportGraph *G_s;
    static int *edge_marks;
    static int *best_tour_nodes;
  };

  friend class PSEP_SimpleDP;
  std::vector<std::list<std::unique_ptr<SimpleTooth> > > light_teeth;
  std::vector<std::list<std::unique_ptr<SimpleTooth> > > heavy_teeth;

  std::vector<int> &best_tour_nodes;
  std::vector<int> &edge_marks;

};

#endif
