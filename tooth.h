#ifndef PSEP_TOOTH_H
#define PSEP_TOOTH_H

#include "Graph.h"

class PSEP_CandTooth {
 public:
 PSEP_CandTooth(std::vector<int> & _tour_nodes, SupportGraph & _G) :
  best_tour_nodes(_tour_nodes), G_s(_G) {
    SimpleTooth::ncount = _tour_nodes.size();
    SimpleTooth::G_s = &_G;
  }

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
    }

    int body_size();
    bool body_contains(const int node_perm);
    static bool C_body_subset(const SimpleTooth &T, const SimpleTooth &R);
    void complement();
    void increment_slack(double *lhs_p, double *rhs_p);


  
  private:
    int root;
    int body_start;
    int body_end;
    double slack;
    bool sandwich;

    int node_index;

    friend class PSEP_CandTooth;
    static int ncount;
    static SupportGraph *G_s;
  };

  std::vector<int> &best_tour_nodes;
  SupportGraph &G_s;

};

#endif