/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/** @file
 * @brief FAST BLOSSOM SEPARATION HEURISTICS
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef PSEP_FASTBLOSSOM_H
#define PSEP_FASTBLOSSOM_H

#include "Graph.hpp"
#include "cc_lpcuts.hpp"
#include "datagroups.hpp"

namespace PSEP {
namespace Cut {

class FastBlossoms {
public:
  FastBlossoms(PSEP::Data::GraphGroup &_graph_dat,
	       PSEP::Data::BestGroup &_best_dat,
	       PSEP::Data::SupportGroup &_supp_dat,
	       PSEP::TourGraph &_TG,
	       PSEP::Cut::LPcutIn &_cutq);

  bool find_cuts();

protected:
  PSEP::Data::GraphGroup &graph_dat;
  PSEP::Data::BestGroup &best_dat;
  PSEP::Data::SupportGroup &supp_dat;
  
  PSEP::TourGraph &TG;
  
  PSEP::Cut::LPcutIn &cutq;

private:
  bool find_odd();
  bool find_gh();
};

}
}

#endif
