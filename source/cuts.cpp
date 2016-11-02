#include "cuts.hpp"

#include <map>

using namespace std;

namespace PSEP {
  
template<>
void CutQueue<HyperGraph>::push_front(const HyperGraph &H)
{
  cut_q.push_front(H);
  if(cut_q.size() > q_capacity){
    cut_q.back().delete_refs();
    cut_q.pop_back();
  }
}

template<>
void CutQueue<HyperGraph>::push_back(const HyperGraph &H)
{
  if(cut_q.size() >= q_capacity){
    cut_q.back().delete_refs();
    cut_q.pop_back();
  }
  cut_q.push_back(H);
}

template<>
void CutQueue<HyperGraph>::pop_front()
{
  cut_q.front().delete_refs();
  cut_q.pop_front();
}

}

using namespace PSEP;

int CutTranslate::get_sparse_row(const HyperGraph &H, vector<int> &rmatind,
				 vector<double> &rmatval, char &sense,
				 double &rhs)
{
  int rval = 0;
  sense = 'G';

  map<int, double> coef_map;
  vector<int> body_nodes;
  int deltacount = 0;
  rmatind.clear();
  rmatval.clear();

  
  rval = H.source_setbank->extract_nodelist(*(H.set_refs[0]),
					    body_nodes);
  PSEP_CHECK_RVAL(rval, "Couldn't extract nodelist, ");
  

  
  GraphUtils::get_delta(body_nodes, edges, &deltacount, delta, edge_marks);
  rval = (deltacount == 0);
  PSEP_CHECK_RVAL(rval, "Body nodes gave empty delta, ");
  

  try {
    rmatind.resize(deltacount);
    rmatval.resize(deltacount);
    for(int i = 0; i < deltacount; i++)
      coef_map[delta[i]] = 1.0;
  }
  catch(const std::bad_alloc &){
    rval = 1; PSEP_GOTO_CLEANUP("Out of memory for sparse row, ");
  }

  switch(H.cut_type){
  case HyperGraph::CutType::Segment:    
    rhs = 2;
    break;

  case HyperGraph::CutType::Blossom:
    int num_teeth = H.set_refs.size() - 1;
    rhs = 1 - num_teeth;
    
    //skip the first ref, it is body above
    for(int i = 1; i < H.set_refs.size(); i++){
      vector<int> edge_tooth;
      int edge_index;
      IntPairMap::iterator find_it;
      
      rval = H.source_setbank->extract_nodelist(*(H.set_refs[i]), edge_tooth);
      PSEP_CHECK_RVAL(rval, "Couldn't extract blossom tooth, ");

      if(edge_tooth.size() != 2){
	rval = 1; PSEP_GOTO_CLEANUP("Blossom tooth has " << edge_tooth.size()
				    << "nodes! ");
      }

      
      find_it = edge_lookup.find(IntPair(fmin(edge_tooth[0], edge_tooth[1]),
					 fmax(edge_tooth[0], edge_tooth[1])));
      rval = (find_it == edge_lookup.end());
      PSEP_CHECK_RVAL(rval, "Couldn't find tooth in edge lookup, ");
            
      edge_index = find_it->second;
      coef_map[edge_index] = -1.0;
    }
    
    break;
  }

  {//scoped temporary variable
    int i = 0;
    for(map<int, double>::iterator it = coef_map.begin(); it != coef_map.end();
	it++){
      rmatind[i] = it->first;
      rmatval[i] = it->second;
      i++;
    }
  }

 CLEANUP:
  if(rval){
    cerr << "CutTranslate<HyperGraph>::get_sparse_row failed, row is invalid\n";
    rmatind.clear();
    rmatval.clear();
  }
  return rval;
}

int CutTranslate::get_sparse_row_if(bool &violated, const HyperGraph &H,
				    const vector<double> &x,
				    vector<int> &rmatind,
				    vector<double> &rmatval, char &sense,
				    double &rhs)
{
  int rval = 0;
  violated = false;
  double activity;

  rval = get_sparse_row(H, rmatind, rmatval, sense, rhs);
  PSEP_CHECK_RVAL(rval, "Couldn't get sparse row, ");

  get_activity(activity, x, rmatind, rmatval);

  switch(sense){
  case 'G':
    violated = (activity < rhs);
    break;
  case 'L':
    violated = (activity > rhs);
    break;
  default:
    rval = 1;
    PSEP_GOTO_CLEANUP("Uncaught row sense " << sense << ", ");
  }

 CLEANUP:
  if(rval)
    std::cerr << "CutTranslate<HyperGraph>::get_sparse_row_if failed, "
	      << "row is invalid.\n";

  if(rval || !violated){
    rmatind.clear();
    rmatval.clear();
  }
  
  return rval;
}

int CutTranslate::is_cut_violated(bool &violated, const HyperGraph &H,
				  vector<double> &x)
{
  int rval = 0;
  vector<int> rmatind;
  vector<double> rmatval;
  char sense; double rhs;

  violated = false;

  rval = get_sparse_row_if(violated, H, x, rmatind, rmatval, sense, rhs);
  PSEP_CHECK_RVAL(rval, "Couldn't test violation, ");

 CLEANUP:
  if(rval)
    cerr << "CutTranslate::is_cut_violated failed\n";
  return rval;
}
