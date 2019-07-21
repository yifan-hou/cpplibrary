#ifndef _TIMER_LINUX_H_
#define _TIMER_LINUX_H_

#include <cstdlib>
#include <chrono>

namespace RUT {

class Timer
{
public:
	Timer();
	~Timer();

	void tic();
	double toc(); // return ms

private:
	std::chrono::high_resolution_clock::time_point _t1;
	std::chrono::high_resolution_clock::time_point _t2;
};

}

#endif // _TIMER_LINUX_H_
