#ifndef PSEP_TIMER_HPP
#define PSEP_TIMER_HPP

#include <chrono>
#include <string>

#include <ctime>

namespace PSEP {

/** A class for recording CPU and wall clock time. */
class Timer {
public:
  Timer(); /**< Construct unnamed timer. */
  
  Timer(const std::string &tname); /**< Construct a named timer. */

  /** Construct a named timer with ratio timer. */
  Timer(const std::string &tname,  const Timer *_ratio_timer);
  
  void start();
  void stop();
  void resume();

  /** Reports the elapsed times and ratios if applicable. */
  void report();

private:
  std::string timer_name;
  
  std::chrono::duration<double> wall_elapsed;
  double cpu_elapsed;

  std::chrono::time_point<std::chrono::system_clock> wall_start;
  std::chrono::time_point<std::chrono::system_clock> wall_end;
  
  double cpu_start;
  double cpu_end;

  /** For use by @ref report to report elapsed times as ratio of other timer.
   * If \p ratio_timer is not `nullptr`, then a call to @ref report will
   * also report wall_elapsed as a ratio of `ratio_timer->wall_elapsed`.
   */
  const Timer *ratio_timer;
};

}

#endif
