/*
 * Various helper functions using ffmpeg
 *
 * Copyright 2016 Yahoo Inc.
 * Licensed under the terms of the Apache 2.0 License.
 * See LICENSE file in the project root for terms.
 *
 * Developer: Yale Song (yalesong@yahoo-inc.com)
 */


#ifndef HECATE_FFMPEG_HPP
#define HECATE_FFMPEG_HPP

#include <stdio.h>
#include <string.h>

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "hecate/file_helper.hpp"

#define BUF_S 256
#define BUF_M 512
#define BUF_L 1024

using namespace std;

namespace hecate {
  const string _ffmpeg = hecate::which("ffmpeg");
  const string _ffmarg = "-loglevel quiet -y";
  const int _lfade = 8; // fade-in/out duration (unit: frame)
  
  // Crop video segment & resize
  static inline void ffmpeg_video_crop(string in_file,
                                       string out_file,
                                       string start_time,
                                       string duration,
                                       int out_width_px)
  {
    in_file = escape_space(in_file);
    
    char filter[BUF_S];
    sprintf( filter, "'scale=%d:trunc(ow/a/2)*2'", out_width_px);
    
    char cmd[BUF_L];
    sprintf( cmd, "%s -ss %s -t %s -i %s -strict -2 -vf %s %s %s",
            _ffmpeg.c_str(), start_time.c_str(), duration.c_str(),
            in_file.c_str(), filter, out_file.c_str(), _ffmarg.c_str());
    
    system( cmd );
  };
  
  static inline void ffmpeg_video_concat(string in_filelist,
                                         string out_file)
  {
    char cmd[BUF_L];
    sprintf( cmd, "%s -f concat -i %s -c copy %s %s",
            _ffmpeg.c_str(), in_filelist.c_str(), out_file.c_str(),
            _ffmarg.c_str() );
    system( cmd );
  };
  
  // Audio fade in/out
  static inline void ffmpeg_audio_fade(string in_file,
                                       string out_file,
                                       double video_duration_sec,
                                       double video_fps)
  {
    in_file = escape_space(in_file);
    
    const double afade_sec = (double)2*_lfade/video_fps;
    const double afade_msec = (int)10000*(afade_sec-floor(afade_sec));
    
    double afos  = video_duration_sec - afade_sec; // audio fade-out start
    int afos_ss  = (int) afos; // NOTE: we don't compute modulo this time
    int afos_mss = (int) 10000*(afos - floor(afos));
    
    char filter[BUF_S];
    sprintf( filter, "'afade=t=in:ss=0:d=0.%04d,"
                      "afade=t=out:st=%d.%04d:d=0.%04d'",
            (int)afade_msec, afos_ss, afos_mss, (int)afade_msec);

    char cmd[BUF_L];
    sprintf( cmd, "%s -i %s -af %s %s %s", _ffmpeg.c_str(),
            in_file.c_str(), filter, out_file.c_str(), _ffmarg.c_str() );
    system( cmd );

  };
  
  // Video fade in/out
  static inline void ffmpeg_video_fade(string in_file,
                                       string out_file,
                                       int video_duration,
                                       bool out_only=false)
  {
    in_file = escape_space(in_file);
    
    char filter[BUF_S];
    if( out_only ) {
      sprintf( filter, "'fade=out:%d:%d'",
              video_duration-_lfade, _lfade);
    }
    else {
      sprintf( filter, "'fade=in:0:%d,fade=out:%d:%d'",
              _lfade, video_duration-_lfade, _lfade);
    }
    
    char cmd[BUF_L];
    sprintf( cmd, "%s -i %s -vf %s %s %s", _ffmpeg.c_str(),
            in_file.c_str(), filter, out_file.c_str(), _ffmarg.c_str());
    system( cmd );
  };
  
  
  // Based on http://blog.pkh.me/p/21-high-quality-gif-with-ffmpeg.html
  static inline void ffmpeg_video2gif(string in_file,
                                      string out_file,
                                      string start_time,
                                      string duration,
                                      int out_fps,
                                      int out_width_px)
  {
    in_file = escape_space(in_file);
    
    string out_dir = hecate::get_dir( std::string(out_file) );
    
    char filter[BUF_S];
    sprintf( filter, "fps=%d,scale=%d:-1:flags=lanczos", out_fps, out_width_px );
    
    char palette[BUF_M];
    sprintf( palette, "%s/palatte.png", out_dir.c_str() );
    
    char time_setup[BUF_S] = "";
    if( !start_time.empty() && !duration.empty() ) {
      sprintf( time_setup, "-ss %s -t %s",
              start_time.c_str(), duration.c_str() );
    }

    char cmd[BUF_L];
    
    // Create a palatte
    sprintf( cmd, "%s %s -i %s -vf '%s,palettegen=stats_mode=diff' %s %s",
            _ffmpeg.c_str(), time_setup,  in_file.c_str(), filter,
            _ffmarg.c_str(), palette );
    system( cmd );
    
    // Convert segment to gif
    sprintf( cmd, "%s %s -i %s -i %s -lavfi '%s [x]; [x][1:v] paletteuse' %s %s",
            _ffmpeg.c_str(), time_setup, in_file.c_str(), palette, filter,
            _ffmarg.c_str(), out_file.c_str());
    system( cmd );
    
    // Delete palette
    sprintf( cmd, "rm %s", palette );
    system( cmd );
  };
  
  
  
}
#endif

