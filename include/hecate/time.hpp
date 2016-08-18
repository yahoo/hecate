/*
 * Time util
 *
 * Copyright 2016 Yahoo Inc.
 * Licensed under the terms of the Apache 2.0 License.
 * See LICENSE file in the project root for terms.
 *
 * Developer: Yale Song (yalesong@yahoo-inc.com)
 */

#ifndef HECATE_TIME_HPP
#define HECATE_TIME_HPP

#include <stdio.h>

#include <ctime>      // display date and time
#include <cmath>      // floor
#include <chrono>
#include <string>

namespace hecate {
  
  typedef std::chrono::high_resolution_clock Clock;
  typedef std::chrono::milliseconds milliseconds;
  
  inline double elapsed_time_ms( Clock::time_point from)
  {
    Clock::time_point now = Clock::now();
    milliseconds ms = std::chrono::duration_cast<milliseconds>( now - from );
    return (double)ms.count();
  }
  
  inline double print_elapsed_time( Clock::time_point from, const char* prefix)
  {
    int msec = (int) elapsed_time_ms(from);
    int sec = msec/1000;
    printf("%s: Elapsed time %02d:%02d:%04d\n", prefix, sec/60, sec%60, msec%1000);
    
    return msec / 1000.0;
  }
  
  inline std::string second2string(double sec, const std::string& format) {
    char buf[128];
    int hh, mm, ss, mss;
    hh  = (int) floor( sec/3600 );
    mm  = (int) floor( sec/60 ) % 60;
    ss  = (int) sec % 60;
    mss = (int) 10000*(sec - floor(sec));
    if( format == "hh:mm:ss.mss" ) {
      sprintf( buf, "%02d:%02d:%02d.%04d", hh,mm,ss,mss );
    }
    else if( format == "mm:ss.mss" ) {
      sprintf( buf, "%02d:%02d.%04d", 60*hh+mm,ss,mss );
    }
    std::string ret = buf;
    return ret;
  }
}
#endif

