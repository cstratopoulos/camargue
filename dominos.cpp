#include<iostream>

#include "dominos.hpp"


using namespace std;
using namespace PSEP;

int Cut<domino>::cutcall(){
  int rval = 0;

  rval = separate();
  if(rval) goto CLEANUP;

  rval = parse_coeffs();
  if(rval) goto CLEANUP;

  rval = add_cut();

 CLEANUP:
  if(rval == 1)
    cerr << "Problem in Cut<domino>::cutcall\n";
  return rval;
}

int Cut<domino>::separate(){
  int rval = 0;

  rval = SimpleDP.separate();
  if(rval) goto CLEANUP;

 CLEANUP:
  if(rval == 1)
    cerr << "Problem in Cut<domino>::separate\n";
  return rval;
}

int Cut<domino>::parse_coeffs(){
  int rval = 0;
  agg_coeffs.resize(PSEPlp_numcols(&m_lp), 0);
  rhs = 0;
  double lhs = 0;
  int RHS;
  int ncount = SimpleDP.light_nodes.size();
  int ecount = SimpleDP.cut_ecap.size();
  vector<int> cut_node_marks(ncount, 0);
  vector<int> domino_delta(ecount, 0);
  int deltacount = 0;

  GraphUtils::get_delta(SimpleDP.cut_nodes.size(),
		     &(SimpleDP.cut_nodes)[0], ecount,
		     &(SimpleDP.cut_elist)[0], &deltacount,
		     &domino_delta[0], &cut_node_marks[0]);

  SimpleDP.parse_domino(deltacount, domino_delta, agg_coeffs, &rhs);

  rhs /= 2;
  RHS = (int) rhs;
  rhs = RHS;

  for(int j = 0; j < agg_coeffs.size(); j++){
    if(((int) agg_coeffs[j]) % 2 == 1){
      rval = 2;
      goto CLEANUP;
    }

    agg_coeffs[j] /= 2;
  }

  for(int j = 0; j < m_lp_edges.size(); j++)
    lhs += m_lp_edges[j] * agg_coeffs[j];

  if(lhs <= rhs){
    rval = 2;
    goto CLEANUP;
  } else{
    cout << "Found simpleDP with lhs/rhs: " << lhs << "/" << rhs << "\n";
  }

 CLEANUP:
  if(rval == 2)
    cerr << "Cut<domino>::parse_coeffs found a bad inequality\n";
  return rval;
}

int Cut<domino>::add_cut(){
  int rval = 0, newrows = 1, newnz;
  vector<int> rmatind;
  vector<double> rmatval;
  char sense[1];
  int rmatbeg[1];

  rmatbeg[0] = 0;
  sense[0] = 'L';

  for(int i = 0; i < agg_coeffs.size(); i++){
    if(agg_coeffs[i] != 0.0){
      rmatind.push_back(i);
      rmatval.push_back(agg_coeffs[i]);
    }
  }
  newnz = rmatind.size();

  rval = PSEPlp_addrows(&m_lp, newrows, newnz, &rhs, sense, rmatbeg,
			&rmatind[0], &rmatval[0]);

  if(rval)
    cerr << "Problem in Cut<domino>::add_cut\n";
  agg_coeffs.clear();
  rhs = 0;

  return rval;
}
