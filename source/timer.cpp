#include "timer.hpp"

#include <iostream>
#include <iomanip>
#include <stdexcept>

#include <cstdio>

using std::cout;
using std::cerr;
using std::setw;
using std::setprecision;
using std::string;
using std::chrono::system_clock;

using std::runtime_error;
using std::exception;

namespace CMR {

Timer::Timer() try
    : timer_name("Unnamed timer"),  wall_elapsed(0), cpu_elapsed(0),
      ratio_timer(nullptr)
{} catch (const exception &e) {
    cerr << e.what() << "\n";
    throw runtime_error("Timer constructor failed.");
}

Timer::Timer(const string &tname) try :
  timer_name(tname), wall_elapsed(0), cpu_elapsed(0), ratio_timer(nullptr) {}
catch (const exception &e) {
    cerr << e.what() << "\n";
    throw runtime_error("Timer constructor failed.");
}

Timer::Timer(const string &tname, const Timer *_ratio_timer) try :
  timer_name(tname), wall_elapsed(0), cpu_elapsed(0), ratio_timer(_ratio_timer)
{} catch (const exception &e) {
    cerr << e.what() << "\n";
    throw runtime_error("Timer constructor failed.");
}

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

void Timer::report(bool show_cpu) const
{
    string name_pad;

    if (timer_name.length() <= 20)
        name_pad = string(20 - timer_name.length(), ' ');
    cout << name_pad << timer_name << ": ";

    printf("%.6fs wall \t", wall_elapsed.count());
    if (ratio_timer)
        printf("(%04.1f%% of %s)",
               (100 * (wall_elapsed.count() /
                       ratio_timer->wall_elapsed.count())),
               ratio_timer->timer_name.c_str());
    printf("\n");

    if (show_cpu) {
        string cpu_pad(22, ' ');
        cout << cpu_pad;
        printf("%.6fs CPU \t", cpu_elapsed);
        if (ratio_timer)
            printf("(%04.1f%% of %s)",
                   (100 * (cpu_elapsed / ratio_timer->cpu_elapsed)),
                   ratio_timer->timer_name.c_str());
        printf("\n");
    }
}

}
