#include "tooth.h"

int PSEP_CandTooth::SimpleTooth::ncount;
SupportGraph * PSEP_CandTooth::SimpleTooth::G_s;


int PSEP_CandTooth::SimpleTooth::body_size(){
  int add_one = (int) (!sandwich); //total is + 1 if tooth is not sandwiched

  return (body_start <= body_end) ? (body_end - body_start + add_one) :
    ((ncount - body_start) + body_end + add_one);
}

bool PSEP_CandTooth::SimpleTooth::body_contains(const int perm_node){
  if(perm_node == root)
    return false;

  return (body_start <= body_end) ?
    (body_start <= perm_node && perm_node <= body_end) :
    (body_start <= perm_node || perm_node <= body_end);
}
