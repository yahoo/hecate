/*
 * HECTATE Video Processing Library
 *
 * Copyright 2016 Yahoo Inc.
 * Licensed under the terms of the Apache 2.0 License.
 * See LICENSE file in the project root for terms.
 *
 * Developer: Yale Song (yalesong@yahoo-inc.com)
 */

#ifndef HECATE_HPP
#define HECATE_HPP

#include <stdio.h>
#include <stdlib.h>
#include <dirent.h>
#include <getopt.h>
#include <sys/stat.h> // mkdir
#include <unistd.h>   // access() function

#include <ctime>      // display date and time
#include <iostream>   // for standard I/O
#include <fstream>    // for file I/O
#include <string>     // for strings
#include <iomanip>    // for controlling float print precision
#include <sstream>    // string to number conversion
#include <chrono>
#include <limits>
#include <numeric>

// OpenMP library
#if defined(_OPENMP)
#include <omp.h>
#endif

// Hecate headers
#include "hecate/sort.hpp"
#include "hecate/time.hpp"
#include "hecate/gapstat.hpp"
#include "hecate/knapsack.hpp"
#include "hecate/cc_parser.hpp"
#include "hecate/file_helper.hpp"
#include "hecate/video_parser.hpp"
#include "hecate/ffmpeg_helper.hpp"

// OpenCV library
#include <opencv2/objdetect/objdetect.hpp>
#include <opencv2/ml/ml.hpp>

using namespace std;

// program options
struct hecate_params {
  string in_video;
  string out_dir;
  string caption;
  int step_sz;          // frame subsampling step size
  int njpg;             // number of thumbnail images
  int ngif;             // number of GIFs
  int lmov;             // length of video summary (in seconds)
  int gif_fps;          // gif play speed
  int jpg_width_px;     // thumbnail image width
  int gif_width_px;     // animated GIF width
  int mov_width_px;     // summary video width
  int max_duration;     // maximum length of video to process (in seconds)
  double fltr_begin_sec;// always filter out x-second frames at the beginning
  double fltr_end_sec;  // always filter out x-second frames at the end
  double invalid_wnd;   // window for dropping neighbor frames of low-quality ones
  bool jpg;             // generate thumbnail jpg
  bool gif;             // generate summary gif
  bool mov;             // generate summary movie
  bool gifsum;          // if enabled, combine individual GIFs into one summary
  bool gifall;          // generate all possible gifs (for debugging purpose)
  bool info_shot;       // print shot boundary info
  bool info_keyfrm;     // print key frame indices
  bool prefer_dynamic;  // if enabled, prefer dynamic scene in highlight
  bool gfl;             // run group-fused lasso as part of shot segmentation
  bool fade;            // disable fade-in/out during shot transition
  bool debug;
  bool display;
  
  hecate_params():
  out_dir("./output"),
  caption(""),
  step_sz(1),
  njpg(5),
  ngif(5),
  lmov(15),
  gif_fps(8),
  jpg_width_px(360),
  gif_width_px(360),
  mov_width_px(360),
  max_duration(-1),
  fltr_begin_sec(-1.0),
  fltr_end_sec(-1.0),
  invalid_wnd(0.15),
  jpg(false),
  gif(false),
  mov(false),
  gifsum(false),
  gifall(false),
  info_shot(false),
  info_keyfrm(false),
  prefer_dynamic(true),
  gfl(false),
  fade(false),
  debug(false),
  display(false)
  {};
};

inline void hecate_parse_params(int argc, char** argv, hecate_params& opt)
{
  static struct option long_options[] = {
    {"in_video",        required_argument, 0, 'i'},
    {"out_dir",         required_argument, 0, 'o'},
    {"step",            required_argument, 0, 's'},
    {"njpg",            required_argument, 0, 'n'},
    {"ngif",            required_argument, 0, 'q'},
    {"lmov",            required_argument, 0, 'l'},
    {"gif_fps",         required_argument, 0, 'f'},
    {"jpg_width_px",    required_argument, 0, 'u'},
    {"gif_width_px",    required_argument, 0, 'v'},
    {"mov_width_px",    required_argument, 0, 'w'},
    {"max_duration",    required_argument, 0, 'd'},
    {"fltr_begin_sec",  required_argument, 0, 'a'},
    {"fltr_end_sec",    required_argument, 0, 'b'},
    {"invalid_wnd",     required_argument, 0, 'k'},
    {"generate_jpg",      no_argument, 0, 'J'},
    {"generate_gif",      no_argument, 0, 'G'},
    {"generate_mov",      no_argument, 0, 'M'},
    {"optimize_gif",      no_argument, 0, 'O'},
    {"generate_gifsum",   no_argument, 0, 'S'},
    {"generate_gifall",   no_argument, 0, 'A'},
    {"print_shot_info",   no_argument, 0, 'T'},
    {"print_keyfrm_info", no_argument, 0, 'K'},
    {"prefer_dynamic",    no_argument, 0, 'V'},
    {"gfl",               no_argument, 0, 'B'},
    {"fade",              no_argument, 0, 'F'},
    {"debug",             no_argument, 0, 'D'},
    {"display",           no_argument, 0, 'C'},
    {0,0,0,0}
  };
  
  while( true ) {
    int opt_idx=0;
    int c = getopt_long( argc, argv,
                        "0:1:2:3:4:5:6:7:8:9:"
                        "a:b:c:d:e:f:g:h:i:j:k:l:m:n:o:p:q:r:s:t:u:v:w:x:y:z:"
                        "A:B:C:D:E:F:G:H:I:J:K:L:M:N:O:P:Q:R:S:T:U:V:W:X:Y:Z",
                        long_options, &opt_idx);
    if( c==-1 ) break;
    switch( c ) {
      case 'i': opt.in_video         = optarg; break;
      case 'o': opt.out_dir          = optarg; break;
      case 's': opt.step_sz          = atoi(optarg); break;
      case 'n': opt.njpg             = atoi(optarg); break;
      case 'q': opt.ngif             = atoi(optarg); break;
      case 'l': opt.lmov             = atoi(optarg); break;
      case 'f': opt.gif_fps          = atoi(optarg); break;
      case 'u': opt.jpg_width_px     = atoi(optarg); break;
      case 'v': opt.gif_width_px     = atoi(optarg); break;
      case 'w': opt.mov_width_px     = atoi(optarg); break;
      case 'd': opt.max_duration     = atof(optarg); break;
      case 'a': opt.fltr_begin_sec   = atof(optarg); break;
      case 'b': opt.fltr_end_sec     = atof(optarg); break;
      case 'k': opt.invalid_wnd      = atof(optarg); break;
      case 'J': opt.jpg              = true; break;
      case 'G': opt.gif              = true; break;
      case 'M': opt.mov              = true; break;
      case 'S': opt.gifsum = opt.gif = true; break;
      case 'A': opt.gifall = opt.gif = true; break;
      case 'T': opt.info_shot        = true; break;
      case 'K': opt.info_keyfrm      = true; break;
      case 'V': opt.prefer_dynamic   = true; break;
      case 'B': opt.gfl              = true; break;
      case 'F': opt.fade             = true; break;
      case 'D': opt.debug            = true; break;
      case 'C': opt.display          = true; break;
    }
  }
  
  // Create output dir (silently fails if dir already exists)
  mkdir( opt.out_dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH );
}



inline void hecate_copyright()
{
  printf("\n");
  printf("-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n");
  printf(" HECATE Copyright 2016 Yahoo Inc.\n");
  printf("   Licensed under the terms of the Apache 2.0 License.\n");
  printf("   Developed by : Yale Song (yalesong@yahoo-inc.com)\n");
  printf("   Built on  : %s %s\n", __TIME__, __DATE__ );
  printf("-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n");
}

inline void hecate_usage()
{
  hecate_params opt;
  printf("USAGE: hecate -i infile [options]\n");
  printf("\n");
  printf("  -i  --in_video      (string)    Input video file\n");
  printf("  -o  --out_dir       (string)    Output directory (%s)\n", opt.out_dir.c_str());
  printf("  -s  --step          (int)       Frame subsampling step size (%d)\n", opt.step_sz);
  printf("  -n  --njpg          (int)       Number of thumbnails to be generated (%d)\n", opt.njpg);
  printf("  -q  --ngif          (int)       Number of GIFs to be generated (%d)\n", opt.ngif);
  printf("  -r  --lmov          (int)       Length of video summary to be generated (in seconds) (%d)\n", opt.lmov);
  printf("  -u  --jpg_width_px  (int)       Pixel width of thumbnail images (%d)\n", opt.jpg_width_px);
  printf("  -v  --gif_width_px  (int)       Pixel width of animated GIFs (%d)\n", opt.gif_width_px);
  printf("  -w  --mov_width_px  (int)       Pixel width of summary video (%d)\n", opt.mov_width_px);
  printf("  --generate_jpg                  Generate thumbnail images\n");
  printf("  --generate_gif                  Generate animated GIFs\n");
  printf("  --generate_mov                  Generate a summary video\n");
  printf("  --generate_gifsum               Generate animated GIFs summary\n");
  printf("  --generate_gifall               Generate all possible animated GIFs\n");
  printf("  --print_shot_info               Print valid shot ranges\n");
  printf("  --print_keyfrm_info             Print keyframe indices\n");
  
  exit(-1);
}



/*******************************************************************************
 
 HECTATE INTERFACE
 
 INPUT:
 hecate_params opt : input option paramaters
 
 OUTPUT:
 vector<int> v_thumb_idx : thumbnail frame indices
 vector<hecate::Range> v_gif_range : highlight shot ranges for GIF creation
 vector<hecate::Range> v_mov_range : highlight shot ranges for video summarization
 
 *******************************************************************************/
void run_hecate( hecate_params& opt,
                 vector<int>& v_thumb_idx,
                 vector<hecate::Range>& v_gif_range,
                 vector<hecate::Range>& v_mov_range );

inline void run_hecate( hecate_params& opt, vector<int>& v_thumb_idx) {
  vector<hecate::Range> v_gif_range, v_mov_range;
  run_hecate( opt, v_thumb_idx, v_gif_range, v_mov_range );
}

inline void run_hecate( hecate_params& opt, vector<hecate::Range>& v_range) {
  vector<int> v_thumb_idx;
  vector<hecate::Range> v_xxx_range;
  
  if( opt.gif )
    run_hecate( opt, v_thumb_idx, v_range, v_xxx_range );
  else if( opt.mov )
    run_hecate( opt, v_thumb_idx, v_xxx_range, v_range );
}


/*******************************************************************************
 
 THUMBNAIL IMAGE GENERATION MODULE
 
 INPUT:
 hecate_params opt        : input option paramaters
 video_metadata meta      : video metadata
 vector<hecate::ShotRanges> v_shot_range
                          : shot boundaries
 const Mat& X             : input features
 const Mat& diff          : n-by-1 vector of frame-by-frame difference scores
 
 OUTPUT:
 vector<int> v_thumb_idx  : vector of frame index numbers for thumbnails
 
 *******************************************************************************/

void detect_thumbnail_frames( hecate_params& opt,
                              hecate::video_metadata& meta,
                              const vector<hecate::ShotRange>& v_shot_range,
                              const Mat& X,
                              const Mat& diff,
                              vector<int>& v_thumb_idx);

void generate_thumbnails( hecate_params& opt, vector<int>& v_thumb_idx );




/*******************************************************************************
 
 VIDEO SUMMARIZATION MODULE
 
 INPUT:
 hecate_params opt        : input option paramaters
 video_metadata meta      : video metadata
 vector<hecate::ShotRanges> v_shot_range
                          : shot boundaries
 const Mat& X             : input features
 const Mat& diff          : n-by-1 vector of frame-by-frame difference scores
 
 OUTPUT:
 vector<hecate::Ranges> v_highlight_range
                          : vector of highlight shot ranges
 
 *******************************************************************************/

void detect_highlight_shots( hecate_params& opt,
                             hecate::video_metadata& meta,
                             const vector<hecate::ShotRange>& v_shot_range,
                             const Mat& X,
                             const Mat& diff,
                             vector<hecate::Range>& v_highlight_range );

void generate_highlight_clips( hecate_params& opt,
                               vector<hecate::Range>& v_highlight_range );


////////////////////////////////////////////////////////////////////////////////
//
// VARIOUS HELPER FUNCTIONS
//
////////////////////////////////////////////////////////////////////////////////


inline void mark_invalid( vector<bool>& vec, int idx, int wnd_sz=0 )
{
  int vec_len = (int)vec.size();
  for(int i=max(0,idx-wnd_sz); i<=min(vec_len-1,idx+wnd_sz); i++)
    vec[i] = false;
}

inline void expand_invalid_frms( vector<bool>& valid, int k )
{
  vector<bool> valid_new = valid;
  
  int nfrm = (int) valid.size();
  for( int pos=1; pos<nfrm-1; pos++ )
  {
    if( !valid[pos] )
      for(int i=max(0,pos-k); i<=min(nfrm-1,pos+k); i++)
        valid_new[i] = false;
  }
  
  valid = valid_new;
}


#endif

