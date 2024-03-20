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

	void tic();
	double toc(); // return ms

private:
	TimePoint _t1;
	TimePoint _t2;
};

}

#endif // _TIMER_LINUX_H_
