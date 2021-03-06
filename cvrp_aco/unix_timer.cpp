/*

       AAAA    CCCC   OOOO   TTTTTT   SSSSS  PPPPP
      AA  AA  CC     OO  OO    TT    SS      PP  PP
      AAAAAA  CC     OO  OO    TT     SSSS   PPPPP
      AA  AA  CC     OO  OO    TT        SS  PP
      AA  AA   CCCC   OOOO     TT    SSSSS   PP

######################################################
##########    ACO algorithms for the TSP    ##########
######################################################

      Version: 1.0
      File:    unix_timer.c
      Author:  Thomas Stuetzle
      Purpose: routines for measuring elapsed time (CPU or real)
      Check:   README.txt and legal.txt
*/

#include <stdio.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <time.h>

#include "timer.h"


static double virtual_time, real_time;
static char format_time[26];


void start_timers(void)
/*    
      FUNCTION:       virtual and real time of day are computed and stored to 
                      allow at later time the computation of the elapsed time 
		      (virtual or real) 
      INPUT:          none
      OUTPUT:         none
      (SIDE)EFFECTS:  virtual and real time are computed   
*/
{
    struct rusage res;
    struct timeval tp;
    
    getrusage( RUSAGE_SELF, &res );
    virtual_time = (double) res.ru_utime.tv_sec +
		   (double) res.ru_stime.tv_sec +
		   (double) res.ru_utime.tv_usec / 1000000.0 +
		   (double) res.ru_stime.tv_usec / 1000000.0;

    gettimeofday( &tp, NULL );
    real_time =    (double) tp.tv_sec +
		   (double) tp.tv_usec / 1000000.0;
}



double elapsed_time(TIMER_TYPE type)
/*    
      FUNCTION:       return the time used in seconds (virtual or real, depending on type) 
      INPUT:          TIMER_TYPE (virtual or real time)
      OUTPUT:         seconds since last call to start_timers (virtual or real)
      (SIDE)EFFECTS:  none
*/
{
    struct rusage res;
    struct timeval tp;
    
    if (type == REAL) {
        gettimeofday( &tp, NULL );
        return( (double) tp.tv_sec +
		(double) tp.tv_usec / 1000000.0
		- real_time );
    }
    else {
        getrusage( RUSAGE_SELF, &res );
        return( (double) res.ru_utime.tv_sec +
		(double) res.ru_stime.tv_sec +
		(double) res.ru_utime.tv_usec / 1000000.0 +
		(double) res.ru_stime.tv_usec / 1000000.0
		- virtual_time );
    }

}


char* get_format_time(void)
{
    time_t timer;
    struct tm* tm_info;
    
    time(&timer);
    tm_info = localtime(&timer);
    
    strftime(format_time, 26, "%Y-%m-%d %H:%M:%S", tm_info);
    return format_time;
}



