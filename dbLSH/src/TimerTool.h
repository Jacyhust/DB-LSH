//This file will implement a cross-platform high resolution timer
//by GS
#ifndef _TIMER_TOOL_H_
#define _TIMER_TOOL_H_

#pragma once



//#ifdef WIN32
////under WIN32, we use QueryPerformanceCounter()
//#include <windows.h>
//#include <ctime>
//#else
//#include <windows.h>
//#include <ctime>
//#endif

#if defined(unix) || defined(__unix__)
//under POSIX system, we use clock_gettime()
//remember we have to use linker option "-lrt"
#include <time.h>
#else
#include <windows.h>
#include <ctime>
#endif

using namespace std;

class HighPrecsionTimer
{
public:
	static HighPrecsionTimer *instance;
	static double getTime()
	{

#ifdef WIN32
		LONGLONG clock_count;
		QueryPerformanceCounter((LARGE_INTEGER *)&clock_count);
		return (double)clock_count/getInstance()->clock_freq;
#endif

#if defined(unix) || defined(__unix__)
		timespec clock;
		clock_gettime(CLOCK_MONOTONIC, &clock);
		return (double)clock.tv_sec+(double)clock.tv_nsec*1e-9;
#endif

	}
protected:

#ifdef WIN32
	LONGLONG clock_freq;
#endif

	static HighPrecsionTimer *getInstance()
	{
		if(instance==NULL) instance=new HighPrecsionTimer();		
		return instance;
	}
private:
	HighPrecsionTimer()
	{

#ifdef WIN32
		QueryPerformanceFrequency((LARGE_INTEGER *)&clock_freq);
#endif

	}

	~HighPrecsionTimer()
	{
		if(instance!=NULL) delete instance;
	}
};

//finally our get timer interface
inline double gettime()
{
	return HighPrecsionTimer::getTime();
}

#endif