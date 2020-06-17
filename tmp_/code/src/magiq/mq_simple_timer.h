/* 
 * File:   simple_timer.h
 * Author: fuad
 *
 * Created on November 26, 2017, 9:32 AM
 */

#ifndef SIMPLE_TIMER_H
#define	SIMPLE_TIMER_H

#include <sys/time.h>  


struct timer {
	timeval i_start_, i_stop_;
	void start() {
		gettimeofday(&i_start_, NULL);
	}
	
	void stop() {
		gettimeofday(&i_stop_, NULL);
	}
	
	double interval() {
		double t1 = i_start_.tv_sec + i_start_.tv_usec/1e6;
		double t2 = i_stop_.tv_sec + i_stop_.tv_usec/1e6;
		return (t2-t1);
	}
};


#endif	/* SIMPLE_TIMER_H */

