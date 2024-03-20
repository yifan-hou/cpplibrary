#include "RobotUtilities/TimerLinux.h"

namespace RUT {
Timer::Timer()
{
    _t1 = Clock::now();
    _t2 = Clock::now();
}


Timer::~Timer()
{
}

void Timer::tic()
{
	_t1 = Clock::now();
}

double Timer::toc()
{
	_t2 = Clock::now();
    return double(std::chrono::duration_cast<std::chrono::nanoseconds>(_t2 - _t1).count())/1e6; // milli second
}

}

