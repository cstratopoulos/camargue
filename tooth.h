#ifndef PSEP_TOOTH_H
#define PSEP_TOOTH_H

#include "Graph.h"

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

  void complement();
  bool body_contains(const int node_index);

  void increment_slack(double *lhs_p, double *rhs_p, const SupportGraph &G);

  static bool C_body_subset(const SimpleTooth &T, const SimpleTooth &R);
  
 private:
  int root;
  int body_start;
  int body_end;
  double slack;
  bool sandwich;

  int node_index;

  static int ncount;
};

#endif
