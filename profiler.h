#ifndef PSEP_TIMER_H
#define PSEP_TIMER_H

#include <iostream>

class PSEP_Timer {
 public:
  static void profile_lp_speed();
  
  static void profile_sep_speed();

  static void profile_seg_speed();
  static void profile_2match_speed();
  static void profile_simple_dp_speed();

  static void profile_lp_cuts();

 private:
  static double pivot_time;
  static double opt_time;
  static double copybase_time;
  
  static double match_time;

  
  static double seg_time;

  
  static double dp_time;
  static double dp_lemma1_time;
  static double dp_lemma2_time;
  static double dp_slack_l1_time;
  static double dp_slack_l2_time;
  static double dp_tree_time;
  static double dp_webedge_time;
  static double dp_gomoryhu_time;
  static double dp_tooth_parse_time;
  static double dp_handle_parse_time;
  static double dp_cutgraph_time;

  static int seg_cutcount;
  
  static int blossom_cutcount;
  static int num_bad2_in_e0;
  static int num_bad2_in_e1;
  
  static int dp_cutcount;
  static int dp_bad_count;

  friend class TSP_Solver;
};

#endif
