#include "RobotUtilities/TimerLinux.h"

#include <iostream>
#include <thread>

namespace RUT {
Timer::Timer() {
  _t1 = Clock::now();
  _t2 = Clock::now();
  _next_loop_start_t = Clock::now();
  _loop_duration_s = std::chrono::duration<double>{-1.};
}

Timer::~Timer() {}

TimePoint Timer::now() { return Clock::now(); }

void Timer::tic() { _t1 = Clock::now(); }

double Timer::toc_ms() {
  _t2 = Clock::now();
  return double(std::chrono::duration_cast<std::chrono::nanoseconds>(_t2 - _t1)
                    .count()) /
         1e6;  // milli second
}

bool Timer::set_loop_rate_hz(double hz) {
  if (hz <= 0) {
    std::cerr << "[Timer] Error: rate must be a positive number." << std::endl;
    return false;
  }
  _loop_duration_s = std::chrono::duration<double>{1. / hz};
  std::cout << "[Timer] debug duration set: " << _loop_duration_s.count()
            << std::endl;

  return true;
}

bool Timer::start_timed_loop() {
  _next_loop_start_t = Clock::now();
  return true;
}

bool Timer::sleep_till_next() {
  _next_loop_start_t +=
      std::chrono::duration_cast<std::chrono::nanoseconds>(_loop_duration_s);
  std::this_thread::sleep_until(_next_loop_start_t);
  return true;
}

}  // namespace RUT
