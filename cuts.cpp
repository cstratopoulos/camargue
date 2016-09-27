#include "cuts.hpp"

#include <map>

using namespace std;
using namespace PSEP;
  
template<>
void PSEP::CutQueue<HyperGraph>::push_front(const HyperGraph &H)
{ 
  cut_q.push_front(H);
  if(cut_q.size() > q_capacity){
    cut_q.back().delete_refs();
    cut_q.pop_back();
  }
}

template<>
void PSEP::CutQueue<HyperGraph>::pop_front()
{
  cut_q.front().delete_refs();
  cut_q.pop_front();
}

