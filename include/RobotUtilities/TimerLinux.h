#ifndef _TIMER_LINUX_H_
#define _TIMER_LINUX_H_

#include <cstdlib>
#include <chrono>

namespace RUT {

typedef std::chrono::high_resolution_clock Clock;
typedef Clock::time_point TimePoint;

class Timer
{
public:
	Timer();
	~Timer();

	// time measurement functions
	void tic();
	double toc_ms(); // return ms

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

}

#endif // _TIMER_LINUX_H_
