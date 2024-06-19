#ifndef _TIMER_LINUX_H_
#define _TIMER_LINUX_H_

#include <chrono>
#include <cstdlib>

namespace RUT {

typedef std::chrono::steady_clock Clock;
typedef Clock::time_point TimePoint;

class Timer {
 public:
  Timer();
  ~Timer();

  // return epoch time
  TimePoint now();

  // duration measurement
  void tic();
  void tic(const TimePoint& time_point);
  double toc_ms();  // return ms
  TimePoint toc_time_point();  // return TimePoint

  // timed loop
  bool set_loop_rate_hz(double hz);
  bool start_timed_loop();
  bool sleep_till_next();

 private:
  TimePoint _t1;
  TimePoint _t2;
  TimePoint _next_loop_start_t;
  std::chrono::duration<double> _loop_duration_s;
};

}  // namespace RUT

#endif  // _TIMER_LINUX_H_
