#ifndef PSEP_TIMER_HPP
#define PSEP_TIMER_HPP

#include <chrono>

#include <ctime>

namespace PSEP {

/** A class for recording CPU and wall clock time. */
class Timer {
public:
  Timer();
  
  void start();
  void stop();
  void resume();

private:
  std::chrono::duration<double> wall_elapsed;
  std::chrono::duration<double> CPU_elapsed;

  std::chrono::time_point<std::chrono::system_clock> wall_start;
  std::chrono::time_point<std::chrono::system_clock> wall_end;
  
  std::clock_t cpu_start;
  std::clock_t cpu_end;
};

}

#endif
