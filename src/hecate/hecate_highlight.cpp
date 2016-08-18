/*
 * HECTATE Video Processing Library - Highlight
 *
 * Copyright 2016 Yahoo Inc.
 * Licensed under the terms of the Apache 2.0 License.
 * See LICENSE file in the project root for terms.
 *
 * Developer: Yale Song (yalesong@yahoo-inc.com)
 */

#include "hecate/hecate.hpp"

using namespace std;
using namespace cv;


void detect_highlight_shots( hecate_params& opt, hecate::video_metadata& meta,
                            const vector<hecate::ShotRange>& v_shot_range,
                            const Mat& X, const Mat& diff,
                            vector<hecate::Range>& v_highlight_range)
{
  if( opt.gif && opt.mov ) {
    fprintf( stderr, "Fatal Error in detect_highlight_shots():"
            "opt.gif and opt.mov cannot be true at the same time.\n" );
    exit(-1);
  }
  
  v_highlight_range.clear();

  const int minK = 5;
  const int maxK = 30;

  ////////////////////////////////////////////////////////////////////////////
  //
  // Discard shots that are too short or too static
  //
  ////////////////////////////////////////////////////////////////////////////

  int min_shot_len, max_shot_len, min_num_shot;
  if( opt.mov ) {
    double addrate_min = 0.5 * sqrt(max(0.0,(double)opt.lmov-15.0)/45.0);
    double addrate_max = sqrt(max(0.0,(double)opt.lmov-15.0)/45.0);
    min_shot_len = floor((3.0+addrate_min) * meta.fps / opt.step_sz);
    max_shot_len = round((1.0+addrate_max) * min_shot_len);
    min_num_shot = ceil( opt.lmov / 3.0 );
    if( opt.debug ) {
      printf("detect_highlight_shots(): "
             "min_shot_len=%d, max_shot_len=%d, min_num_shot=%d, "
             "addrate_min=%f, addrate_max=%f\n",
             min_shot_len, max_shot_len, min_num_shot,
             addrate_min, addrate_max);
    }
  }
  else {
    min_shot_len = floor(1.5 * meta.fps / opt.step_sz);
    max_shot_len = ceil (4.5 * meta.fps / opt.step_sz);
    min_num_shot = opt.ngif;
  }

  // Active set containing candidates for highlight shots
  vector<hecate::Range> v_candidates;
  for(size_t i=0; i<v_shot_range.size(); i++) {
    hecate::Range shot(v_shot_range[i].start, v_shot_range[i].end);
    shot.v_idx = v_shot_range[i].v_idx;
    int shot_len = shot.length();
    double avg_diff=0;
    for(int i=shot.start; i<=shot.end; i++)
      avg_diff += diff.at<double>(i);
    avg_diff/=shot.length();
    
    // Discard if too short.
    if( shot_len < min_shot_len ) {
      //printf("Discard SHORT (shot_len=%d, avg_diff=%f) ", shot_len, avg_diff); shot.print();
      continue;
    }
    
    // Discard if too static
    if( avg_diff < 0.05 ) {
      //printf("Discard STATIC (shot_len=%d, avg_diff=%f) ", shot_len, avg_diff); shot.print();
      continue;
    }
    
    // Add the shot
    v_candidates.push_back( shot );
  }

  // If there's not enough shots, merged adjacent shots and add them
  if( (int)v_candidates.size() < min_num_shot ) {
    if( opt.debug ) {
      printf("detect_highlight_shots(): "
             "Not enough candidate shots (current=%d/%d, min=%d)\n",
             (int)v_candidates.size(), (int)v_shot_range.size(), min_num_shot);
    }
    for(size_t i=0; i<v_shot_range.size(); i++) {
      // Merge shots until it's longer than min_shot_len
      size_t j=i;
      int cumlen=0;
      while( (j+1)<v_shot_range.size() ) {
        cumlen += v_shot_range[j].length();
        if( cumlen<min_shot_len ) { j++; } else { break; }
      }
      // Construct a shot
      hecate::Range shot(v_shot_range[i].start, v_shot_range[j].end);
      shot.v_idx.clear();
      for(size_t k=i; k<=j; k++)
        shot.v_idx.insert( shot.v_idx.end(),
            v_shot_range[k].v_idx.cbegin(), v_shot_range[k].v_idx.cend());

      // Check if this shot is already in v_candidates
      bool exists = false;
      for( size_t k=0; k<v_candidates.size(); k++ ) {
        if( v_candidates[k].start >= shot.start &&
           v_candidates[k].start <= shot.end ) {
          exists = true; break;
        }
      }
      if( !exists && shot.length() >= min_shot_len ) {
        v_candidates.push_back( shot );
      }
      i = j+1;
    }

    // Sort candidate shots in a chronogical order
    vector<hecate::Range> v_candidates_tmp;
    vector<int> v_shot_nsrt; // shots in non-sorted order
    for( size_t i=0; i<v_candidates.size(); i++ )
      v_shot_nsrt.push_back( v_candidates[i].start );
    // sort wrt frame index in an ascending order
    vector<size_t> v_srt_idx; // contains idx
    vector<int> v_srt_val; // contains value (start frame idx)
    hecate::sort( v_shot_nsrt, v_srt_val, v_srt_idx );
    for(size_t i=0; i<v_srt_idx.size(); i++) {
      v_candidates_tmp.push_back( v_candidates[v_srt_idx[i]] );
    }
    v_candidates = v_candidates_tmp;
  }
  
  
  // Shorten shots if too long.
  // Prefer static shot; potentially reduces jitter
  for(size_t i=0; i<v_candidates.size(); i++) {
    hecate::Range shot = v_candidates[i];
    int shot_len = shot.length();
    if( shot_len > max_shot_len ) {
      while( shot_len > max_shot_len ) {
        Scalar mu, sigma, mu2, sigma2;
        meanStdDev( diff(Rect(0,shot.start+1,1,shot.length()-1)), mu, sigma );
        meanStdDev( diff(Rect(0,shot.start,  1,shot.length()-1)), mu2, sigma2 );
        if( sigma[0] < sigma2[0] )
          shot.start++;
        else
          shot.end--;
        shot_len = shot.length();
      }
    }
    v_candidates[i] = shot;
  }
  
  
  ////////////////////////////////////////////////////////////////////////////
  //
  // Group visually simimlar shots
  //
  // Represent each shot as a piece-wise constant multi-dimensional signal
  // Then run kmeans (k=#nshots_highlight).
  //
  ////////////////////////////////////////////////////////////////////////////
  
  // Prepare data for kmeans
  // Compute piece-wise constant feature representation
  Mat km_data( (int)v_candidates.size(), X.cols, X.type() );
  for( size_t shotid=0; shotid<v_candidates.size(); shotid++ ) {
    hecate::Range s = v_candidates[shotid];
    Mat Xsub = X( Rect(0, s.start, X.cols, s.length()) );
    cv::reduce( Xsub, km_data.row(shotid), 0, CV_REDUCE_AVG );
  }
  
  // Perform k-means (repeat 5 times)
  Mat km_lbl; // integer row vector; stores cluster IDs for every sample.
  Mat km_ctr; // one row per each cluster center.
  int km_k = min(maxK, min((int)v_candidates.size(), max(minK, min_num_shot)));
  hecate::perform_kmeans( km_data, km_lbl, km_ctr, km_k, 5 );
  
  // measure cluster size
  vector<int> v_shotlen;
  for(size_t i=0; i<v_candidates.size(); i++) {
    v_shotlen.push_back( v_candidates[i].length() );
  }
  vector<int> clust_sz(km_k,0);
  for(int i=0; i<km_lbl.rows; i++)
    clust_sz[ km_lbl.at<int>(i) ] += v_shotlen[i];
  
  // sort wrt cluster size in an ascending order
  vector<size_t> v_srt_idx; // contains cluster id
  vector<int> v_srt_val;    // contains cluster size
  hecate::sort( clust_sz, v_srt_val, v_srt_idx );
  
  
  ////////////////////////////////////////////////////////////////////////////
  //
  // Shot evaluation
  //
  // Criteria: prefer the most dynamic shot among others within the same cluster
  //
  ////////////////////////////////////////////////////////////////////////////
  
  // Pre-compute per-shot diff score avg & stddev
  Mat shot_diff_avg(v_candidates.size(), 1, diff.type(), Scalar(0,0,0)[0]);
  Mat shot_diff_std(v_candidates.size(), 1, diff.type(), Scalar(0,0,0)[0]);
  for( size_t shotid=0; shotid<v_candidates.size(); shotid++ ) {
    Scalar mu, sigma;
    hecate::Range s = v_candidates[shotid];
    meanStdDev( diff(Rect(0,s.start,1,s.length())), mu, sigma );
    shot_diff_avg.at<double>(shotid) = (double) mu[0];
    shot_diff_std.at<double>(shotid) = (double) sigma[0];
  }

  vector<int> v_shot_len;
  vector<double> v_shot_score;
  for(size_t i=0; i<v_candidates.size(); i++) {
    v_shot_len.push_back( v_candidates[i].length() );
    v_shot_score.push_back( 0.0 );
  }

  // Start with the largest cluster
  for(int i=0; i<km_k; i++)
  {
    int k = v_srt_idx[km_k-i-1];
    hecate::Range best_range(-1,-1);
    double best_val = -numeric_limits<double>::max();
    
    // Collect IDs & scores of the shots in the k-th cluster
    vector<int> v_tmp_shotid;
    vector<double> v_tmp_score;
    for(int shotid=0; shotid<km_lbl.rows; shotid++) {
      int lbl = km_lbl.at<int>(shotid);
      if( lbl==k ) {
        double val = shot_diff_avg.at<double>(shotid);
        if( val > best_val ) {
          best_range = v_candidates[shotid];
          best_val = val;
        }
        v_tmp_shotid.push_back( shotid );
        v_tmp_score.push_back( val );
      }
    }
    
    // If GIF mode, add shots here
    if( opt.gif ) {
      if( best_range.start>=0 && best_range.length()>0 )
        v_highlight_range.push_back( best_range );
      if( (int)v_highlight_range.size() >= opt.ngif )
        break;
    }
    
    // If MOV mode, store shot scores
    if( opt.mov ) {
      // Sort wrt avg_diff score, ascending order
      vector<size_t> v_srt_idx2;
      vector<double> v_srt_val2;
      hecate::sort( v_tmp_score, v_srt_val2, v_srt_idx2 );
      
      // Record normalized scores
      for( size_t i=0; i<v_srt_idx2.size(); i++ ) {
        // if prefer dynamic, higher avg_diff gets priority
        int order = (opt.prefer_dynamic) ? v_srt_idx2.size()-1-i : i;
        v_shot_score[v_tmp_shotid[v_srt_idx2[order]]] = 100.0/(1.0+(double(i)));
      }
    }
  }
  
  
  ////////////////////////////////////////////////////////////////////////////
  //
  // Shot selection
  //
  // Criteria: prefer the most dynamic shot among others within the same cluster
  //
  ////////////////////////////////////////////////////////////////////////////
  
  // Return all available segments
  if ( opt.gif ) {
    if( opt.gifall ) {
      v_highlight_range = v_candidates;
    }
  }
  
  if( opt.mov ) {
    int budget = floor(opt.lmov * meta.fps / opt.step_sz);
    budget += (int)(0.5 * max_shot_len); // buffer
    
    // solve 0/1 knapsack
    vector<bool> sol;
    hecate::solve_01knapsack( v_shot_score, v_shot_len, budget, sol );
    
    // merge shots if only they are only half a second apart
    int merge_thrsh = ceil(0.5 * meta.fps / opt.step_sz);
    for( size_t i=0; i<sol.size(); i++ ) {
      if( sol[i] ) {
        hecate::Range r = v_candidates[i];
        if( v_highlight_range.empty() ) {
          v_highlight_range.push_back( v_candidates[i] );
        }
        else {
          int dist = r.start - v_highlight_range[v_highlight_range.size()-1].end + 1;
          if( dist < merge_thrsh )
            v_highlight_range[v_highlight_range.size()-1].end = r.end;
          else
            v_highlight_range.push_back( v_candidates[i] );
        }
      }
    }
  }

  for(size_t i=0; i<v_highlight_range.size(); i++) {
    v_highlight_range[i].start *= opt.step_sz;
    v_highlight_range[i].end *= opt.step_sz;
  }
  
  
  ////////////////////////////////////////////////////////////////////////////
  //
  // Post Processing
  //
  // Ensures summary video length is exactly the same as required
  //
  ////////////////////////////////////////////////////////////////////////////
  
  // Ensures highlight ranges sums up to exactly target length
  if( opt.mov ) {
    int target = round(opt.lmov * meta.fps);
    int curlen = 0;
    for(size_t i=0; i<v_highlight_range.size(); i++) {
      curlen += v_highlight_range[i].length();
    }
    //printf("target=%d, cur=%d\n", target, curlen);
    if( curlen!=target ) {
      // If summary is too long
      while( curlen>target ) {
        // Pick a random shot
        int idx = rand() % (int) v_highlight_range.size();
        hecate::Range shot = v_highlight_range[idx];
        // Remove frame that reduces stddev of the shot by the most
        Scalar mu, sigma, mu2, sigma2;
        meanStdDev( diff(Rect(0,shot.start+1,1,shot.length()-1)), mu, sigma );
        meanStdDev( diff(Rect(0,shot.start,1,shot.length()-1)), mu2, sigma2 );
        if( sigma[0] < sigma2[0] )
          v_highlight_range[idx].start++;
        else
          v_highlight_range[idx].end--;
        curlen--;
      }
      // If video is longer than target, but summary is too short
      int shot_edge_buf = 10;
      int failed_attempt = 0;
      while( curlen<target ) {
        // Pick a random shot
        int idx = rand() % (int) v_highlight_range.size();
        hecate::Range shot = v_highlight_range[idx];
        
        // Add frame that increases stddev of the shot by the least
        bool failed = false;
        if( v_highlight_range[idx].start<shot_edge_buf ) {
          if( idx+1<v_highlight_range.size() && v_highlight_range[idx+1].start >
             v_highlight_range[idx].end+1+shot_edge_buf ) {
              v_highlight_range[idx].end++;
          }
          else {
            failed = true;
          }
        }
        else if( v_highlight_range[idx].end>meta.nframes-shot_edge_buf ) {
          if( idx>0 && v_highlight_range[idx-1].end <
             v_highlight_range[idx].start-1-shot_edge_buf ) {
            v_highlight_range[idx].start--;
          }
          else {
            failed = true;
          }
        }
        else {
          Scalar mu, sigma, mu2, sigma2;
          meanStdDev( diff(Rect(0,shot.start-1,1,shot.length()+1)), mu, sigma );
          meanStdDev( diff(Rect(0,shot.start,1,shot.length()+1)), mu2, sigma2 );
          if( sigma[0] < sigma2[0] ) {
            if( v_highlight_range[idx-1].end <
               v_highlight_range[idx].start-1-shot_edge_buf ) {
              v_highlight_range[idx].start--;
            }
            else {
              failed = true;
            }
          }
          else {
            if( v_highlight_range[idx+1].start >
               v_highlight_range[idx].end+1+shot_edge_buf) {
              v_highlight_range[idx].end++;
            }
            else {
              failed = true;
            }
          }
        }
        if( !failed )
          curlen++;
        else
          failed_attempt++;
        if( failed_attempt>10*(int)v_highlight_range.size() ) {
          if( opt.debug ) {
            printf("Failed to meet the target length in video summary\n");
          }
          break;
        }
      }
    }
    //curlen = 0;
    //for(size_t i=0; i<v_highlight_range.size(); i++)
    //  curlen += v_highlight_range[i].length();
    //printf("target=%d, cur=%d\n", target, curlen);
  }
}


////////////////////////////////////////////////////////////////////////////
//
// Highlight clip generation
//
////////////////////////////////////////////////////////////////////////////

void generate_highlight_clips( hecate_params& opt, vector<hecate::Range>& v_highlight_range )
{
  if( opt.gif && opt.mov ) {
    fprintf( stderr, "Fatal Error in generate_highlight_clips():"
            "opt.gif and opt.mov cannot be true at the same time.\n" );
    exit(-1);
  }
  
  // prefix for hiddden files generated during execution
  const char *cdot = "__tmp__";
  
  string filename = hecate::get_filename( std::string(opt.in_video) );
  
  VideoCapture vr( opt.in_video );
  double fps = vr.get(CV_CAP_PROP_FPS);
  vr.release();
  
  // Sort shots in chronological order
  vector<int> v_shot_nsrt; // shots in non-sorted order
  for( size_t i=0; i<v_highlight_range.size(); i++ )
    v_shot_nsrt.push_back( v_highlight_range[i].start );
  
  // sort wrt frame index in an ascending order
  vector<size_t> v_srt_idx; // contains idx
  vector<int> v_srt_val; // contains value (start frame idx)
  hecate::sort( v_shot_nsrt, v_srt_val, v_srt_idx );
  
  // For concatenating vieos clips
  char filelist[512];
  sprintf( filelist, "%s/%s%s_seg.txt",
          opt.out_dir.c_str(), cdot, filename.c_str() );
  FILE *ptr_filelist = fopen( filelist, "w+" );
  
  // Parse ttml-format caption, if provided
  vector<hecate::tcc> ttml;
  if( !opt.caption.empty() ) {
    parse_ttml( opt.caption, ttml );
    if( ttml.empty() ) {
      fprintf( stderr, "generate_animated_gifs: "
              "Caption file %s cannot be read\n", opt.caption.c_str() );
    }
    // Convert second to frame index
    for(size_t i=0; i<ttml.size(); i++) {
      ttml[i].start *= fps;
      ttml[i].end *= fps;
    }
  }
  
  // GENERATE CLIPS
  char cmd[1024];
  char infile[512];
  char outfile[512];
  
  for( size_t i=0; i<v_highlight_range.size(); i++ )
  {
    int shotid = v_srt_idx[i];
    hecate::Range r = v_highlight_range[shotid];
    double sec_from = (double) (r.start+1) / fps;
    double sec_duration = (double) (r.end-r.start) / fps;
    string start_pos = hecate::second2string( sec_from, "hh:mm:ss.mss" );
    string duration = hecate::second2string( sec_duration, "hh:mm:ss.mss" );
    
    // --------------------- VIDEO CLIP GENERATOR ---------------------
    if( opt.mov )
    {
      // Crop video segment
      sprintf( outfile, "%s/%s%s_seg%03d.mp4",
              opt.out_dir.c_str(), cdot, filename.c_str(), shotid);
      hecate::ffmpeg_video_crop( opt.in_video, std::string(outfile),
                             start_pos, duration, opt.mov_width_px );
      
      // Apply video fade-in/out
      sprintf( infile, "%s/%s%s_seg%03d.mp4",
              opt.out_dir.c_str(), cdot, filename.c_str(), shotid);
      sprintf( outfile, "%s/%s%s_segV%03d.mp4",
              opt.out_dir.c_str(), cdot, filename.c_str(), shotid);
      if( opt.fade ) {
        hecate::ffmpeg_video_fade( std::string(infile), std::string(outfile),
                               r.end-r.start+1,i==0 );
      }
      else {
        sprintf( cmd, "mv %s %s", infile, outfile ); system( cmd );
      }
      
      // Apply audio fade-in/out
      sprintf( infile, "%s/%s%s_segV%03d.mp4",
              opt.out_dir.c_str(), cdot, filename.c_str(), shotid);
      sprintf( outfile, "%s/%s%s_segAV%03d.mp4",
              opt.out_dir.c_str(), cdot, filename.c_str(), shotid);
      hecate::ffmpeg_audio_fade( std::string(infile), std::string(outfile),
                             sec_duration, fps );
      
      // Log filename for concat
      fprintf( ptr_filelist, "file %s%s_segAV%03d.mp4\n",
              cdot, filename.c_str(), shotid);
    }
    // --------------------- VIDEO CLIP GENERATOR ---------------------
    
    // --------------------- GIF CLIP GENERATOR ---------------------
    if( opt.gif ) {
      if( opt.gifsum )
      {
        // Crop video segment
        sprintf( outfile, "%s/%s%s_seg%03d.mp4",
                opt.out_dir.c_str(), cdot, filename.c_str(), shotid);
        hecate::ffmpeg_video_crop( opt.in_video, std::string(outfile),
                               start_pos, duration, opt.mov_width_px );
        
        // Convert video to gif
        sprintf( infile, "%s", outfile );
        sprintf( outfile, "%s/%s_%02d.gif",
                opt.out_dir.c_str(), filename.c_str(), shotid);
        hecate::ffmpeg_video2gif( std::string(infile), std::string(outfile),
                              "", "", opt.gif_fps, opt.gif_width_px );
        
        // Log filename for concat
        fprintf( ptr_filelist, "file %s%s_seg%03d.mp4\n",
                cdot, filename.c_str(), shotid);
      }
      else
      {
        // Convert video to gif, direct memory access without cropping
        sprintf( outfile, "%s/%s_%02d.gif",
                opt.out_dir.c_str(), filename.c_str(), shotid);
        hecate::ffmpeg_video2gif( opt.in_video, std::string(outfile),
                              start_pos, duration, opt.gif_fps,
                              opt.gif_width_px );
      }
    }
    // --------------------- GIF CLIP GENERATOR ---------------------
  }
  
  // Close filelist
  fclose( ptr_filelist );
  
  if( opt.mov )
  {
    // Concatenate segments
    sprintf( outfile, "%s/%s_sum.mp4",
            opt.out_dir.c_str(), filename.c_str() );
    hecate::ffmpeg_video_concat( filelist, outfile );
  }
  
  if( opt.gif && opt.gifsum )
  {
    // Concatenate segments
    sprintf( outfile, "%s/%s%s_segsum.mp4",
            opt.out_dir.c_str(), cdot, filename.c_str() );
    hecate::ffmpeg_video_concat( filelist, outfile );
    
    // Convert video to gif
    sprintf( infile, "%s", outfile );
    sprintf( outfile, "%s/%s_sum.gif",
            opt.out_dir.c_str(), filename.c_str());
    hecate::ffmpeg_video2gif( std::string(infile), std::string(outfile),
                          "", "", opt.gif_fps, opt.gif_width_px );
  }
  
  // Clean up
  sprintf( cmd, "rm %s/%s%s_seg*.*",
          opt.out_dir.c_str(), cdot, filename.c_str() );
  system( cmd );
}

