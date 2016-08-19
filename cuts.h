/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 *                                                                            
 *                CUT TEMPLATE AND CUT STRUCTURE DEFINITIONS                  
 *                                                                            
 * An empty class template for the Cut data type used by the CutControl class 
 * and structure definitions for each type of cut used therein.
 * A specialization for a given cut type should support the named operations, 
 * while including additional functions/data as needed. 
 *                                                                            
 *    CLASS TEMPLATES:                                                        
 *      Cut<cut_type>
 *        public:
 *          int cut_call()
 *          Wrapper function calling private member functions in order
 *        private:
 *          int separate(), parse_coeffs(), add_cut()
 *          separate - the cut-specific (primal) separation routine which
 *                     stores a cut, if found, in the cut-specific storage
 *                     object. See specific headers for details. 
 *                   - returns 0 if successful,
 *                             1 if an error occurs,
 *                             2 if no cut is found
 *          parse_coeffs() - turns cut structure into coefficients that can be
 *                           added to the LP.
 *                         - returns 0 if successful, 1 if an error occurs
 *          add_cut()    - adds the cut found to the LP
 *                       - returns 0 if successful, 1 if an error occurs
 *        
 *    STRUCTURES:
 *      seg
 *        int start, int end, double viol
 *      Stores segment cuts -- subtour cuts arising from segments of the 
 *      incumbent best tour
 *      start, end - the start and endpoints of the segment, indicated as 
 *                   indices for accessing the vector best_tour_nodes. 
 *                   Thus the associated segment is best_tour_nodes[start] to
 *                   best_tour_nodes[end]
 *      viol - the amount by which the cut is violated
 *
 *      

\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef PSEP_CUTS_H
#define PSEP_CUTS_H

#include <memory>
#include <vector>

#include "tooth.h"

namespace PSEP {
  
  struct seg {
  seg(int _start, int _end, double _slack) :
    start(_start), end(_end), viol(_slack) {}
    int start;
    int end;
    double viol;

    bool operator <(const seg &val) const {
      return viol > val.viol;
    }
  };

  struct blossom {
  blossom(std::vector<int> &_handle, int _cut_edge, double _val) :
    handle(_handle), cut_edge(_cut_edge), cut_val(_val){}

    bool operator< (const blossom &val) const {
      return cut_val < val.cut_val;
    }

    std::vector<int> handle;
    int cut_edge;
    double cut_val;
  };

  struct domino {
    domino(){}
    domino(std::vector<int> &_handle,
	   std::vector<std::shared_ptr<CandTooth::SimpleTooth>> _teeth) :
    handle(_handle), used_teeth(_teeth) {}

    std::vector<int> handle;
    std::vector<std::shared_ptr<CandTooth::SimpleTooth>> used_teeth;
  };

  struct general {
  general(bool _gomory, bool _disj, bool _mir) :
    gomory_frac(_gomory), disjunctive(_disj), rounding(_mir) {}

    bool gomory_frac;
    bool disjunctive;
    bool rounding;
  };

  
  template<typename cut_t>
    class Cut {
  public:
    //int cut_call(){ return 1;}

  private:
    /* int separate(){ return 1;} */
    /* int parse_coeffs(){ return 1;} */
    /* int add_cut(){return 1;} */
  };
}

#endif
