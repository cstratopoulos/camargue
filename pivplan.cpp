#include "pivplan.h"

using namespace PSEP;


PivotPlan::PivotPlan(int _ncount, bool _branch,
		     std::vector<ParamPair> ParamList) :
  current_edge_ratio(_ncount),
  ncount(_ncount),
  bash_on_regardless(false),
  branch(_branch),
  max_time(0),
  max_augrounds(0),
  goal_ratio(1) {
  while(!ParamList.empty()){
    switch(ParamList.back().first)
      {
      case Params::TimeLimit:
	max_time = ParamList.back().second;
	std::cout << "TIMELIMIT: set to " << max_time << "\n";
	break;
      case Params::AugLimit:
	max_augrounds = ParamList.back().second;
	std::cout << "AUGLIMIT: set to " << max_augrounds << "\n";
	break;
      case Params::EdgeRatio:
	goal_ratio = ParamList.back().second;
	std::cout << "GOAL RATIO: set to " << goal_ratio << "\n";
	break;
      }	
    ParamList.pop_back();
  }
}

PivotPlan::PivotPlan(int _ncount, Presets preset) :
  current_edge_ratio(_ncount),
  ncount(_ncount),
  bash_on_regardless(false),
  branch(false),
  max_time(0),
  max_augrounds(0),
  goal_ratio(1){
  switch(preset){
  case Presets::BASH_ON:
    bash_on_regardless = true;
    break;
  case Presets::ROOT:
    goal_ratio = 2;
    max_augrounds = 300;
    break;
  case Presets::BRANCH:
    std::cout << "PRESET: BRANCH\n";
    branch = true;
    max_augrounds = 250;
    break;
  case Presets::SPARSE:
    break;
  }
}
