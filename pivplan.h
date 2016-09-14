/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 *
 *                PIVOT PLAN PRESETS AND CLASS DEFINITION
 *
 * This file contains the PivotPlan class, and some enumerations for presets.
 * The PivotPlan allows the PureCut solver to be written in sufficient
 * generality that it may be used on its own, or embedded in the ABC solver.
 * PivotPlan contains info about when to 'give up' on the solution process
 * and whether or not to perform edge elimination.
 *
 * TODO: Add an output interval parameter (freq of print statements)
 *       Maybe figure out variadic initialization list, or interactive init.
 *       Maybe move functions for reporting on pivoting rounds here, to clean
 *       up the PureCut control flow
 *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
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
    /*
     * Presets are preset pivot plans, used most often in practice, as opposed
     * to tuning parameters manually. See Params just below for what the 
     * parameters mean.
     *
     * ROOT - intended to be used to solve the root node of an ABC call
     *    EdgeRatio -
     *    AugLimit -
     * BRANCH - intended to be used to solve a B&B node where some clamping
     *    has been enforced.
     *    EdgeRatio -
     *    AugLimit -
     * BASH_ON - bash on regardless: solving proceeds until either optimality
     *     is proved or cut generation fails
     */
    enum class Presets { ROOT, BRANCH, BASH_ON };

    /*
     * Params are parameters that can be used to define/initialize a PivotPlan
     * They describe conditions under which the solver may prematurely 
     * terminate without proving optimality
     * Zero or negative values indicate that the criterion will not be used
     * See below for functions testing each of these conditions
     *
     * TimeLimit - solver will terminate after this much time elapses
     * AugLimit - solver will terminate after this many rounds without
     *    pivoting to an improved tour
     * EdgeRatio - solver will terminate if enough edges are eliminated so that
     *     the ratio of edges/nodes is within +/- 1 of EdgeRatio
     *     This has no effect unless used in conjunction with a plan 
     *     that performs regular edge elimination
     */
    enum class Params {TimeLimit, AugLimit, EdgeRatio};
    typedef std::pair<Params, int> ParamPair;


    /* The default constructor creates a BASH_ON plan */
  PivotPlan() : current_edge_ratio(INFINITY),
      ncount(1), bash_on_regardless(true), branch(false) {}
    /* Constructor may also accept a Preset or vector of ParamPairs --
     * the parameters to be used and their values. With both of these, we
     * must pass the problem node count to the plan so that edge ratios
     * can be computed. 
     */
    PivotPlan(int _ncount, Presets Preset);
    PivotPlan(int _ncount, bool _branch, std::vector<ParamPair> ParamList);

    
    void start_timer() {start_time = PSEP::zeit();}

    /*
     * condition is used as the termination for the main while() loop in
     * PureCut. Given a Preset or parameters, it will evaluate to true
     * for as long as the associated criteria are satisfied
     */
    bool condition(const int augrounds) {
      if(bash_on_regardless) return true;

      return(UnderTimeLimit() && UnderAugLimit(augrounds) && AboveGoalRatio());
    }

    /*
     * This bool returns true if reduced cost edge fixing/elimination is to
     * be performed under the current PivotPlan
     */
    bool perform_elim(){
      return goal_ratio > 1 || bash_on_regardless;
    }

    /* This function will print the reason why termination occurred */
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

    /* current_edge_ratio is the ratio of edges in the graph to ncount */
    double current_edge_ratio;
    const int ncount;
    /* returns whether or not this is a branching plan */
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
