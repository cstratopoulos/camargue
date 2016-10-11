#include "fastblossoms.hpp"

#include <iostream>

using namespace std;
using namespace PSEP;

int Cut<fastblossom>::separate()
{
  int rval = 0;
  int component_count = 0;
  int *component_sizes = (int *) NULL;
  int *component_nodes = (int *) NULL;
}

int Cut<fastblossom>::add_cuts()
{
  return 1;
}

int Cut<fastblossom>::cutcall()
{
  cout << "The fastblossom cutcall" << endl;
  int rval = separate();
  if(rval) goto CLEANUP;

  rval = add_cuts();
  if(rval) goto CLEANUP;

 CLEANUP:
  if(rval == 1)
    cerr << "Problem in Cuts<fastblossom>::cutcall\n";
  frac_indices.clear();
  frac_elist.clear();
  frac_ecap.clear();
  return rval;
}
