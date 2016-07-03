#include "profiler.h"

double PSEP_Timer::pivot_time = 0;
double PSEP_Timer::opt_time = 0;
double PSEP_Timer::copybase_time = 0;
  
double PSEP_Timer::match_time = 0;
double PSEP_Timer::seg_time = 0;
    
double PSEP_Timer::dp_time = 0;
double PSEP_Timer::dp_lemma1_time = 0;
double PSEP_Timer::dp_lemma2_time = 0;
double PSEP_Timer::dp_slack_l1_time = 0;
double PSEP_Timer::dp_slack_l2_time = 0;
double PSEP_Timer::dp_tree_time = 0;
double PSEP_Timer::dp_webedge_time = 0;
double PSEP_Timer::dp_gomoryhu_time = 0;
double PSEP_Timer::dp_tooth_parse_time = 0;
double PSEP_Timer::dp_handle_parse_time = 0;
double PSEP_Timer::dp_cutgraph_time = 0;

int PSEP_Timer::seg_cutcount = 0;
int PSEP_Timer::blossom_cutcount = 0;
int PSEP_Timer::num_bad2_in_e0 = 0;
int PSEP_Timer::num_bad2_in_e1 = 0;
int PSEP_Timer::dp_cutcount = 0;
int PSEP_Timer::dp_bad_count = 0;
