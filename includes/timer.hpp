#ifndef PSEP_TIMER_HPP
#define PSEP_TIMER_HPP

#include <chrono>
#include <string>

namespace PSEP {

/** A class for recording CPU and wall clock time. */
class Timer {
public:
  Timer(); /**< Construct unnamed timer. */
  
  Timer(const std::string &tname); /**< Construct a named timer. */

  /** Construct a named timer with ratio timer.
   * In the presence of a ratio_timer, calls to @ref profile will report
   * elapsed time(s) as a fraction of those from ratio_timer.
   */
  Timer(const std::string &tname,  const Timer *_ratio_timer);

  /** Start accumulating time at current time, resetting elapsed times. */
  void start();

  /** Stop accumulating times and add them to elapsed times. */
  void stop();

  /** Start accumulating time without resetting elapsed times. */
  void resume();

  /** Reports the elapsed times and ratios if applicable. 
   * If \p show_cpu is true, cpu_elapsed will be reported as well. If not,
   * only wall_elapsed will be reported. Generally show_cpu should be `false`
   * unless the process spawns multiple threads. 
   */
  void report(bool show_cpu);

private:

  /** The name of the Timer.
   * timer_name will be used by @ref report to identify the time(s) printed. 
   * For most nicely formatted output, keep timer_name at 20 chars or less.
   */
  std::string timer_name;

  /** Elapsed wall clock or stopwatch time. */
  std::chrono::duration<double> wall_elapsed;

  /** Elapsed CPU time. */
  double cpu_elapsed;

  /** Wall clock start time. */
  std::chrono::time_point<std::chrono::system_clock> wall_start;

  /** Wall clock end time. */
  std::chrono::time_point<std::chrono::system_clock> wall_end;

  /** CPU clock start time. */
  double cpu_start;

  /** CPU clock end time. */
  double cpu_end;

  /** For use by @ref report to report elapsed times as ratio of other timer.
   * If \p ratio_timer is not `nullptr`, then a call to @ref report will
   * also report wall_elapsed and cpu_elapsed as a ratio of those respective
   * values from ratio_timer. 
   */
  const Timer *ratio_timer;
};

}

#endif
