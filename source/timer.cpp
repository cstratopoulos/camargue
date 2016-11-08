#include "timer.hpp"

#include <iostream>

using std::cout;
using std::string;
using std::chrono::system_clock;

namespace PSEP {

Timer::Timer() :
  wall_elapsed(0), cpu_elapsed(0), ratio_timer(nullptr) {}

Timer::Timer(const string &tname) :
  timer_name(tname), wall_elapsed(0), cpu_elapsed(0), ratio_timer(nullptr) {}

Timer::Timer(const string &tname, const Timer *_ratio_timer):
  timer_name(tname), wall_elapsed(0), cpu_elapsed(0), ratio_timer(_ratio_timer)
{}

void Timer::start()
{
  wall_elapsed = std::chrono::duration<double>(0);
  cpu_elapsed = 0;

  resume();
}

void Timer::stop()
{
  cpu_end = std::clock();
  wall_end = system_clock::now();

  wall_elapsed += wall_end - wall_start;
  cpu_elapsed += (cpu_end - cpu_start) / CLOCKS_PER_SEC;
}

void Timer::resume()
{  
  cpu_start = std::clock();
  wall_start = system_clock::now();
}

void Timer::report()
{

  cout << "    ";
  if(!timer_name.empty())
    cout << timer_name;
  if(ratio_timer)
    cout << " (part of " << ratio_timer->timer_name << ")";
  if(!timer_name.empty());
  cout << ":\n";
  
  cout << "        " << wall_elapsed.count() << "s wall ";
  if(ratio_timer)
    cout << "(ratio "
	 << (wall_elapsed.count() / ratio_timer->wall_elapsed.count())
	 <<  ")";
  cout << "\n";

  cout << "        " << cpu_elapsed << "s CPU ";
  if(ratio_timer)
    cout << "(ratio "
	 << (cpu_elapsed / ratio_timer->cpu_elapsed)
	 << ")";
  cout << "\n";
}

}
