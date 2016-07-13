#ifndef PSEP_SIMPLEDP_H
#define PSEP_SIMPLEDP_H

#include "tooth.h"

class PSEP_SimpleDP {
 public:
  PSEP_SimpleDP(std::vector<int> & _tour_nodes, SupportGraph & _G,
		std::vector<int> & _edge_marks) :
  G_s(_G), candidates(_tour_nodes, _G, _edge_marks) {}
		
  
  void test_build_collection();

 private:
  SupportGraph &G_s;
  PSEP_CandTooth candidates;
};

#endif
