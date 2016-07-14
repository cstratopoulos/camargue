#include "simpleDP.h"
using namespace std;

void PSEP_SimpleDP::test_build_collection(){
  candidates.build_collection();
  int ncount = G_s.node_count;
  int light_total = 0, heavy_total = 0;

  cout << "Printing light teeth...\n";
  for(int i = 0; i < ncount; i++){
    light_total += candidates.light_teeth[i].size();
  }

  cout << light_total << " light teeth in total\n";

  cout << "Printing heavy teeth...\n";
  for(int i = 0; i < ncount; i++){
    heavy_total += candidates.heavy_teeth[i].size();
  }

  cout << heavy_total << " heavy teeth in total\n";

  cout << "Total number of teeth with slack less than one: "
       << light_total + heavy_total << endl;

  cout << "Ratio of light teeth to n^2 is : "
       << ((double) light_total / (ncount * ncount)) << "\n";

  cout << "Ratio of heavy teeth to n^3 is : "
       << ((double) heavy_total / (ncount * ncount * ncount)) << "\n";

  cout << "Ratio of total number of teeth to n^3 is: "
       << ((double) (light_total + heavy_total) / (ncount * ncount * ncount))
       << "\n";

}
