/*
 * Video Parser
 *
 * Copyright 2016 Yahoo Inc.
 * Licensed under the terms of the Apache 2.0 License.
 * See LICENSE file in the project root for terms.
 *
 * Developer: Yale Song (yalesong@yahoo-inc.com)
 */

#ifndef HECATE_VIDEO_PARSER_HPP
#define HECATE_VIDEO_PARSER_HPP

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

// OpenCV library
#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>

#if CV_MAJOR_VERSION > 3
#include <opencv2/videoio/legacy/constants_c.h>
#endif

#include "hecate/sort.hpp"
#include "hecate/gflseg.hpp"
#include "hecate/gapstat.hpp"
#include "hecate/shot_range.hpp"
#include "hecate/hist_opencv.hpp"
#include "hecate/image_metrics.hpp"

using namespace std;
using namespace cv;

namespace hecate {
  
  // video metadata struct
  struct video_metadata {
    int nframes;
    int width;
    int height;
    double duration; // in seconds
    double fps;
  };
  
  struct parser_params {
    int step_sz;
    double fltr_begin_sec;
    double fltr_end_sec;
    double max_duration; // in seconds
    bool gfl;      // use group-fused lasso to refine boundaries
    bool fltr_lq;  // filter low-quality frames
    bool fltr_rdt; // filter redundant frames
    bool debug;
    bool ignore_rest; // if video is too long, ignore the rest (true) or
                      // adjust step_sz (false). for vidtag this should be
                      // false, while for hecate this should be true
    
    parser_params():
    step_sz(1),
    fltr_begin_sec(0),
    fltr_end_sec(0),
    max_duration(-1),
    gfl(false),
    fltr_lq(true),
    fltr_rdt(true),
    debug(false),
    ignore_rest(false)
    {};
  };
  
  inline void print_video_metadata( const string filename,
                                    const hecate::video_metadata m ) {
    printf("%s\n  seconds=%.2f, nframes=%d, fps=%.2f, resolution=%dx%d\n",
           filename.c_str(), m.duration, m.nframes, m.fps, m.width, m.height);
  };
  
  class VideoParser {
  public:
    VideoParser();
    
    /* Perform video parsing */
    // max_duration is there to handle too long videos
    //   (100K frms ~= 1 hr long with 30 fps)
    // filter_first and filter_last is there to filter out
    //   a few frames from the beginning and at the end
    //   to manually handle logos, ending credits, etc.
    vector<hecate::ShotRange> parse_video(const string& in_video,
                                       hecate::parser_params opt);
    
    /* Get the number of valid frames */
    int get_nfrm_valid();
    
    /* Display video with shot boundary information */
    void play_video_filtered(const string& in_video,
                             int step_sz=1,
                             int max_frm_len=360);
    
    /* Get a vector of booleans representing valid frames */
    inline void get_frm_valid(vector<bool>& vec) {vec=_v_frm_valid;}
    
    /* Get effective step size */
    inline int get_effective_step_size() {return _step_sz;}
    
    /* Get extracted features */
    inline const Mat get_frame_features() {return _X_feat;}
    inline const Mat get_frame_diff_features() {return _X_diff;}
    
    
  private:
    /* Read video into the memory */
    int read_video(const string& in_video,
                   int step_sz=1,
                   double max_duration=30*60,
                   bool ignore_rest=false,
                   int max_frm_len=160);
    
    /* Filter first end last few frames */
    void filter_heuristic(double fltr_begin_sec=0,
                          double fltr_end_sec=0);
    
    /* Filter out low-quality frames */
    void filter_low_quality(double thrsh_bright=0.075,
                            double thrsh_sharp=0.08,
                            double thrsh_uniform=0.80 );
    
    /* Filter out frames during transition */
    void filter_transition(double thrsh_diff=0.50,
                           double thrsh_ecr=0.10 );
    
    /* Filter out redundant frames */
    void filter_redundant_and_obtain_subshots();
    
    /* Update shot ranges, filter out shots if too short */
    void update_shot_ranges(int min_shot_len=5);

    /* Perform post processing. Break up shots if too long */
    void post_process(double min_shot_sec=2.0, bool gfl=false);
    
    /* Perform SBD using heuristics: add next big diff if
     the new shot length is longer than min_shot_len */
    void sbd_heuristic(vector<double> v_diff, vector<int>& jump,
                       int njumps, int min_shot_len );
    
    // Helper functions
    
    /* Extract pyramid of histogram features */
    void extract_histo_features(int pyr_level=2,
                                bool omit_filtered=true,
                                int nbin_color_hist=128,
                                int nbin_edge_ori_hist=8,
                                int nbin_edge_mag_hist=8);
    
    void mark_invalid( vector<bool>& vec, int idx, int wnd_sz=0 );
    void mark_invalid( vector<bool>& vec, vector<string>& vec2,
                      int idx, const string msg, int wnd_sz=0 );
    
    void release_memory();
    
    
  public:
    hecate::video_metadata meta;
    vector<bool> _v_frm_valid;    // filtered frames
    vector<string> _v_frm_log;    // filtered frames msgs (debug)
    
  private:
    bool _debug;
    bool _display;
    int _step_sz;
    int _nfrm_total;   // number of frames BEFORE sampling
    int _nfrm_given;   // number of frames AFTER sampling
    int _video_width;
    int _video_height;
    double _video_fps;
    double _video_sec;
    
    Mat _X_feat; // frame-wise feature representation
    Mat _X_diff; // n-by-1 frame-by-frame difference first-order derivative
    Mat _X_ecr;  // n-by-1 ecr first-order derivative
    
    vector<Mat> _v_frm_rgb;       // rgb frames
    vector<Mat> _v_frm_gray;      // gray-scale frames
    
    vector<hecate::ShotRange> _v_shot_ranges;
  };
}


#endif
