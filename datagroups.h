#ifndef PSEP_DATAGROUP_H
#define PSEP_DATAGROUP_H

#include<vector>
#include<memory>

#include<math.h>

#include "PSEP_util.h"
#include "Graph.h"
#include "lp.h"

namespace PSEP {
  namespace Data {
    struct GraphGroup {
      GraphGroup(char *fname, std::unique_ptr<CCdatagroup> &dat);
  
      Graph m_graph;
      std::vector<int> island;
      std::vector<int> delta;
      std::vector<int> edge_marks;
    };

    struct BestGroup {
      BestGroup(const Graph &graph, std::unique_ptr<CCdatagroup> &dat);
      std::vector<int> best_tour_edges;
      std::vector<int> best_tour_nodes;
      std::vector<int> perm;

      double min_tour_value;
    };

    struct LPGroup {
      LPGroup(const Graph &m_graph, PSEP::LP::Prefs &_prefs,
	      const std::vector<int> &perm);
      ~LPGroup(){PSEPlp_free(&m_lp);}
      PSEPlp m_lp;  
      std::vector<double> m_lp_edges;
      std::vector<int> old_colstat;
      std::vector<int> old_rowstat;
      PSEP::LP::Prefs prefs;
    };

    struct SupportGroup {
      SupportGraph G_s;
      std::vector<int> support_indices;
      std::vector<int> support_elist;
      std::vector<double> support_ecap;
    };
  }
}

#endif
