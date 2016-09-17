#include "cutcall.h"

using namespace std;
using namespace PSEP;

int CutControl::primal_sep(const int augrounds, const LP::PivType stat){
  int rval = 0;
  
  int segval = 2, matchval = 2, dpval = 2;
  double segtime, matchtime, dptime;

  segtime = zeit();
  segval = segments.cutcall();
  if(segval == 1){
    rval = 1;
    goto CLEANUP;
  }
  
  //
  // if(!segval) cout << "Added segment cut, row number "
  // 		   << (PSEPlp_numrows(&m_lp) - 1) << "\n";
  //
  
  segtime = zeit() - segtime;
  total_segtime += segtime;
  total_segcalls++;

  matchtime = zeit();
  matchval = blossoms.cutcall();
  if(matchval == 1){
    rval = 1;
    goto CLEANUP;
  }

  //
  // if(!matchval) cout << "Added 2match cut, row number "
  // 		   << (PSEPlp_numrows(&m_lp) - 1) << "\n";
  //
  matchtime = zeit() - matchtime;
  total_2mtime += matchtime;
  total_2mcalls++;

  if(prefs.dp_threshold >= 0 && augrounds > 0 &&
     (augrounds % prefs.dp_threshold == 0)){
    if(segval == 2 && matchval == 2 && stat != LP::PivType::Subtour){
      dptime = zeit();
      dpval = dominos.cutcall();
      if(dpval == 1){
	rval = 1;
	goto CLEANUP;
      }
      dptime = zeit() - dptime;
      total_dptime += dptime;	
    }
  }

  if(segval == 2 && matchval == 2 && dpval == 2)
    rval = 2;

 CLEANUP:
  if(rval == 1)
    cerr << "Problem in CutControl::primal_sep: ";
  if(segval == 1)
    cerr << "Cuts<seg>::cutcall\n";
  if(matchval == 1)
    cerr << "Cuts<blossom>::cutcall\n";
  if(dpval == 1)
    cerr << "Cuts<domino>::cutcall\n";
  return rval;
}

int CutControl::safe_gomory_sep(){
  double gentime = zeit();
  int rval = safe_gomory.cutcall();
  gentime = zeit() - gentime;
  total_gentime += gentime;
  total_gencalls++;
  if(rval == 1)
    cerr << "Problem in CutControl::safe_gomory_sep\n";
  return rval;
}
