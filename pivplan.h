#ifndef PSEP_PIVPLAN_H
#define PSEP_PIVPLAN_H

#include<cmath>
#include<vector>
#include<utility>
#include<iostream>
#include<iomanip>

#include "PSEP_util.h"

namespace PSEP {
  class PivotPlan {
  public:
    enum class Presets { ROOT, BRANCH, BASH_ON };
    enum class Params {TimeLimit, AugLimit, EdgeRatio};
    typedef std::pair<Params, int> ParamPair;


  PivotPlan() : current_edge_ratio(INFINITY),
      ncount(1), bash_on_regardless(true), branch(false) {}
    PivotPlan(int _ncount, Presets Preset);
    PivotPlan(int _ncount, bool _branch, std::vector<ParamPair> ParamList);

    
    void start_timer() {start_time = PSEP::zeit();}
    
    bool condition(const int augrounds) {
      if(bash_on_regardless) return true;

      return(UnderTimeLimit() && UnderAugLimit(augrounds) && AboveGoalRatio());
    }

    bool perform_elim(){
      return goal_ratio > 1 || bash_on_regardless;
    }

    void profile(const int augrounds) {
      if(!UnderTimeLimit()){
	std::cout << "artificial time limit of "
		  << max_time <<"s (runtime: "
		  << runtime << "s)\n";
	  return;
      }

      if(!UnderAugLimit(augrounds)){
	std::cout << "limit of "
		  << max_augrounds << " rounds no augmentation\n";
	return;
      }

      if(!AboveGoalRatio()){
	std::cout << "goal edges/ncount ratio of "
		  << goal_ratio << " (currently: "
		  << std::setprecision(2) << current_edge_ratio
		  << std::setprecision(6) << ")\n";
	return;
      }
    }
    
    double current_edge_ratio;
    const int ncount;
    bool is_branch() {return branch;}

  private:
    bool bash_on_regardless;
    bool branch;

    bool UnderTimeLimit(){
      if(max_time <= 0) return true;
      runtime = PSEP::zeit() - start_time;
      return runtime  < max_time;
    }
    double runtime;
    double start_time;
    double max_time;

    bool UnderAugLimit(const int augrounds){
      if(max_augrounds <= 0) return true;
      return augrounds < max_augrounds;
    }
    int max_augrounds;

    bool AboveGoalRatio(){
      if(goal_ratio <= 1) return true;
      return (current_edge_ratio > goal_ratio &&
	      fabs(current_edge_ratio - goal_ratio) >= 1.01);
    }
    int goal_ratio;
  };

  typedef PivotPlan::Presets PivPresets;
  typedef PivotPlan::Params PivParams;
  typedef std::pair<PivParams, int> PivParamPair;
}

#endif
