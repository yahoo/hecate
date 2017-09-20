/*
 * Video Parser
 *
 * Copyright 2016 Yahoo Inc.
 * Licensed under the terms of the Apache 2.0 License.
 * See LICENSE file in the project root for terms.
 *
 * Developer: Yale Song (yalesong@yahoo-inc.com)
 */

#include "hecate/video_parser.hpp"

using namespace std;
using namespace cv;
using namespace hecate;

/*-----------------------------------------------------------------------*/
VideoParser::VideoParser()
/*-----------------------------------------------------------------------*/
{
  _debug = false;
  _display = false;
  
  _nfrm_total = 0;
  _nfrm_given = 0;
  
  _v_shot_ranges.clear();
}

/*-----------------------------------------------------------------------*/
void VideoParser::release_memory()
/*-----------------------------------------------------------------------*/
{
  _v_frm_rgb.clear();
  _v_frm_gray.clear();
}

/*-----------------------------------------------------------------------*/
vector<hecate::ShotRange> VideoParser::parse_video(const string& in_video,
                                                hecate::parser_params opt)
/*-----------------------------------------------------------------------*/
{
  _debug = opt.debug;
  
  int ret = read_video( in_video, opt.step_sz, opt.max_duration,
                       opt.ignore_rest );
  if( ret<0 ) {
    fprintf( stderr, "VideoParser: Failed to open input video: %s\n",
            in_video.c_str());
    return vector<hecate::ShotRange>();
  }
  
  // Record metadata
  meta.nframes  = _nfrm_total;
  meta.width    = _video_width;
  meta.height   = _video_height;
  meta.fps      = _video_fps;
  meta.duration = _video_sec;
  
  // Frame filtering
  if( opt.fltr_begin_sec>.0 || opt.fltr_end_sec>.0 )
    filter_heuristic(opt.fltr_begin_sec, opt.fltr_end_sec);
  
  if( opt.fltr_lq )
    filter_low_quality();
  
  filter_transition();
  
  // Extract feature representation
  extract_histo_features();
  
  // Post-process (break up shots if too long)
  double min_shot_len_sec = 2.0;
  post_process(min_shot_len_sec, opt.gfl);
  
  release_memory();
  
  // Store shot information
  update_shot_ranges();
  
  // Subshot detection; subshot info is stored during this step
  if( opt.fltr_rdt )
    filter_redundant_and_obtain_subshots();
  
  if( _display )
    play_video_filtered( in_video, 2, 360 );
  
  // Produce result
  return _v_shot_ranges;
}

/*------------------------------------------------------------------------
 Read input video and store as RGB and GRAY formats.
 ------------------------------------------------------------------------*/
int VideoParser::read_video( const string& in_video, int step_sz,
                            double max_duration, bool ignore_rest,
                            int max_frm_len)
/*-----------------------------------------------------------------------*/
{
  if( _debug )
    printf("VideoParser: read_video(\"%s\", \n\tstep_sz=%d, "
           "max_duration=%.2f, max_frm_len=%d, ignore_rest=%d)\n",
           in_video.c_str(), step_sz, max_duration, max_frm_len, ignore_rest);
  
  VideoCapture vr( in_video );
  if( !vr.isOpened() ) {
    return -1;
  }
  
  _nfrm_total   = vr.get(CV_CAP_PROP_FRAME_COUNT);
  _video_width  = vr.get(CV_CAP_PROP_FRAME_WIDTH);
  _video_height = vr.get(CV_CAP_PROP_FRAME_HEIGHT);
  _video_fps    = max(1.0, vr.get(CV_CAP_PROP_FPS));
  if( _video_fps!=_video_fps )
    _video_fps = 29.97;
  
  int max_nfrms = (max_duration<0)
  ? _nfrm_total : round(max_duration*_video_fps);
  
  // need to store _step_sz for computing min_shot_len in sbg_gflseg()
  _step_sz = step_sz;
  if( step_sz>10 ) {
    _step_sz = 10;
    fprintf( stderr, "VideoParser: The maximum step size is 10"
            " (provided is %d)\n", step_sz );
  }
  // if video is too long, and ignore_rest is false, adjust step size
  if( !ignore_rest && max_nfrms > 0 && _nfrm_total > max_nfrms ) {
    _step_sz = ceil( _nfrm_total/max_nfrms );
    fprintf( stderr, "VideoParser: Video too long (%d frames),"
            "increasing step size to %d\n", _nfrm_total, step_sz );
  }
  _step_sz = max(1, _step_sz);
  
  int maxlen = max( _video_width, _video_height );
  double rsz_ratio = (maxlen>max_frm_len)
  ? (double) max_frm_len / maxlen : -1.0;
  
  // Read RGB frames
  _nfrm_total=0; // reset because vr.get is known to be inaccurate
  while( true )
  {
    Mat frm;
    vr >> frm;
    if( frm.empty() ) break;
    
    if( _nfrm_total % _step_sz == 0 ) {
      if( rsz_ratio>0 )
        resize( frm, frm, Size(), rsz_ratio, rsz_ratio, CV_INTER_LINEAR );
      _v_frm_rgb.push_back( frm );
      
      // if video is too long, and ignore_rest is true, cut the rest
      if( ignore_rest && _nfrm_total>=max_nfrms ) {
        if( _debug )
          printf("VideoParser: video too long, "
                 "cutting at (%d)-th frame\n", _nfrm_total);
        break;
      }
    }
    _nfrm_total++;
  }
  vr.release();
  
  _nfrm_given = (int) _v_frm_rgb.size();
  _video_sec  = (double)_nfrm_total/_video_fps;
  
  // RGB2GRAY
  _v_frm_gray.assign( _nfrm_given, Mat() );
  
#pragma omp parallel for
  for( int i=0; i<_nfrm_given; i++ )
  {
    Mat frm_gray;
    cvtColor( _v_frm_rgb[i], frm_gray, CV_BGR2GRAY );
    GaussianBlur( frm_gray, frm_gray, Size(3,3), 0, 0 );
    frm_gray.copyTo( _v_frm_gray[i] );
  }
  
  _v_frm_valid.assign( _nfrm_given, true );
  _v_frm_log.assign( _nfrm_given, " " );
  
  _X_diff = Mat( _nfrm_given, 1, CV_64F, Scalar(0,0,0)[0] );
  _X_ecr  = Mat( _nfrm_given, 1, CV_64F, Scalar(0,0,0)[0] );
  
  return 0;
}


/*------------------------------------------------------------------------
 Manually invalidate first and last frames (heuristic)
 ------------------------------------------------------------------------*/
void VideoParser::filter_heuristic(double fltr_begin_sec, double fltr_end_sec)
/*-----------------------------------------------------------------------*/
{
  int fltr_begin_nfrms = ceil(fltr_begin_sec * _video_fps / (double)_step_sz);
  int fltr_end_nfrms = ceil(fltr_end_sec * _video_fps / (double)_step_sz);
  
  for(int i=0; i<fltr_begin_nfrms; i++)
    mark_invalid(_v_frm_valid, _v_frm_log, i, "[Begin]");
  for(int i=0; i<fltr_end_nfrms; i++)
    mark_invalid(_v_frm_valid, _v_frm_log, _nfrm_total-1-i, "[End]");
}


/*------------------------------------------------------------------------
 Use a collection of heuristics to filter out low-quality frames.
 ------------------------------------------------------------------------*/
void VideoParser::filter_low_quality( double thrsh_bright,
                                     double thrsh_sharp,
                                     double thrsh_uniform )
/*-----------------------------------------------------------------------*/
{
  // filter at most n percent of the total frames
  int nfrm_nperc = (int)(0.15*_nfrm_given);
  
  vector<double> v_brightness(_nfrm_given,0.0);
  vector<double> v_sharpness(_nfrm_given,0.0);
  vector<double> v_uniformity(_nfrm_given,0.0);
  
#pragma omp parallel for
  for( int i=0; i<_nfrm_given; i++ ) {
    v_brightness[i] = hecate::calc_brightness( _v_frm_rgb[i]  );
    v_sharpness[i]  = hecate::calc_sharpness(  _v_frm_gray[i] );
    v_uniformity[i] = hecate::calc_uniformity( _v_frm_gray[i] );
  }
  
  vector<size_t> v_srt_idx; // contains sorted indices
  vector<double> v_srt_val; // contains sorted values
  
  // DARK frame detection
  hecate::sort( v_brightness, v_srt_val, v_srt_idx );
  for( int i=0; i<nfrm_nperc; i++ )
    if( v_srt_val[i] <= thrsh_bright )
      mark_invalid(_v_frm_valid, _v_frm_log, v_srt_idx[i], "[Dark]");
  
  // BLURRY frame detection
  hecate::sort( v_sharpness, v_srt_val, v_srt_idx );
  for( int i=0; i<nfrm_nperc; i++ )
    if( v_srt_val[i] <= thrsh_sharp )
      mark_invalid(_v_frm_valid, _v_frm_log, v_srt_idx[i], "[Blurry]");
  
  // UNIFORM frame detection
  hecate::sort( v_uniformity, v_srt_val, v_srt_idx );
  for( int i=0; i<nfrm_nperc; i++ )
    if( v_srt_val[_nfrm_given-i-1] >= thrsh_uniform )
      mark_invalid(_v_frm_valid, _v_frm_log, v_srt_idx[_nfrm_given-i-1], "[Uniform]");
}


/*------------------------------------------------------------------------
 Frames around shot boundaries are usually low-quality, filter them out
 Use two SBD methods: frame-by-frame difference and ECR
 ------------------------------------------------------------------------*/
void VideoParser::filter_transition( double thrsh_diff, double thrsh_ecr )
/*-----------------------------------------------------------------------*/
{
  int nfrm_nperc = (int)(0.10*_nfrm_given); // n percent of the total frames
  
  vector<double> v_diff(_nfrm_given, 0.0);
  vector<double> v_ecr(_nfrm_given, 0.0);
  
  // sort wrt cluster size in an ascending order
  vector<size_t> v_srt_idx; // contains sorted indices
  vector<double> v_srt_val;    // contains sorted values
  
  // compute the first-order derivative frame-by-frame difference
  int img_sz = _v_frm_gray[0].cols * _v_frm_gray[0].rows;
#pragma omp parallel for
  for( int i=1; i<_nfrm_given-1; i++ ) {
    v_diff[i] = (double) (cv::norm( (_v_frm_rgb[i]-_v_frm_rgb[i-1]) )
                          + cv::norm( (_v_frm_rgb[i]-_v_frm_rgb[i+1]) )) / (2.0*img_sz);
    _X_diff.at<double>(i) = v_diff[i];
  }
  
  // compute edge-change-ratio (ECR)
  {
    int dl_sz = 5; // dilute size
    Mat dl_elm = getStructuringElement(MORPH_CROSS,
                                       Size(2*dl_sz+1, 2*dl_sz+1),
                                       Point(dl_sz, dl_sz));
    
    // Pre-compute edge & edge dilation
    vector<Mat> v_edge( _nfrm_given, Mat() );    // edge images
    vector<Mat> v_edge_dl( _nfrm_given, Mat() ); // edge-diluted images
    
#pragma omp parallel for
    for(int i=0; i<_nfrm_given; i++)
    {
      Mat tmp;
      double theta = threshold(_v_frm_gray[i],tmp,0,255,CV_THRESH_BINARY|CV_THRESH_OTSU);
      Canny( _v_frm_gray[i], v_edge[i], theta, 1.2*theta);
      dilate( v_edge[i], v_edge_dl[i], dl_elm );
      v_edge[i] -= 254; v_edge_dl[i] -= 254;
    }
    
    // Transition detection using ECR (edge change ratio)
#pragma omp parallel for
    for(int i=1; i<_nfrm_given; i++)
    {
      double rho_out, rho_in;
      rho_out = 1.0 - min(1.0,sum(v_edge[i-1].mul(v_edge_dl[i]))[0]/max(1e-6,sum(v_edge[i-1])[0]));
      rho_in  = 1.0 - min(1.0,sum(v_edge_dl[i-1].mul(v_edge[i]))[0]/max(1e-6,sum(v_edge[i-1])[0]));
      
      v_ecr[i] = max(rho_out,rho_in); // edge change ratio
      _X_ecr.at<double>(i) = v_ecr[i];
    }
  }
  
  // CUT detection
  hecate::sort( v_diff, v_srt_val, v_srt_idx );
  for( int i=0; i<nfrm_nperc; i++ )
    if( v_srt_val[_nfrm_given-i-1] >= thrsh_diff )
      mark_invalid(_v_frm_valid, _v_frm_log, v_srt_idx[_nfrm_given-i-1], "[Cut]" );
  
  // TRANSITION detection (cut, fade, dissolve, wipe)
  hecate::sort( v_ecr, v_srt_val, v_srt_idx );
  for( int i=0; i<nfrm_nperc; i++ )
    if( v_srt_val[_nfrm_given-i-1] >= thrsh_ecr )
      mark_invalid(_v_frm_valid, _v_frm_log, v_srt_idx[_nfrm_given-i-1], "[ECR]" );
}



/*------------------------------------------------------------------------
 Perform change point detection either by solving group-fused LASSO of
 "The group fused Lasso for multiple change-point detection" (2011) by
 Kevin Bleakley and Jean-Philippe Vert, or by using heuristics
 -----------------------------------------------------------------------*/
void VideoParser::post_process(double min_shot_sec, bool gfl)
/*-----------------------------------------------------------------------*/
{
  hecate::Segmenter seg;
  int start_idx=-1, end_idx=-1, shotlen=-1;
  int min_shot_len = min_shot_sec * _video_fps / _step_sz;
  int max_shot_len = 3 * min_shot_len;
  double thrsh_gfl = 0.25;
  
  for( size_t i=0; i<_v_frm_valid.size(); i++ )
  {
    if( start_idx<0 && _v_frm_valid[i] ) {
      start_idx = i;
    }
    if( start_idx>=0 && (!_v_frm_valid[i] || i+1==_v_frm_valid.size()) )
    {
      end_idx = i;
      shotlen = end_idx-start_idx+1;
      if( shotlen >= max_shot_len )
      {
        int njumps = floor(shotlen/min_shot_len);
        vector<int> jump;
        
        // Solve group-fused LASSO
        if( gfl ) {
          Mat Xsub( shotlen, _X_feat.cols, _X_feat.type() );
          for( int r=start_idx; r<=end_idx; r++ ) {
            _X_feat.row(r).copyTo( Xsub.row(r-start_idx) );
          }
          Xsub.convertTo(Xsub,CV_64F);
          
          seg.gflseg( Xsub, jump, njumps, thrsh_gfl );
        }
        // Use heuristics
        else {
          vector<double> v_diff;
          for(int i=start_idx; i<=end_idx; i++) {
            v_diff.push_back( _X_diff.at<double>(i) );
          }
          sbd_heuristic( v_diff, jump, njumps, min_shot_len );
        }
        
        for(size_t k=0; k<jump.size(); k++) {
          mark_invalid(_v_frm_valid, _v_frm_log, start_idx+jump[k]-1, "[GFL]" );
          mark_invalid(_v_frm_valid, _v_frm_log, start_idx+jump[k], "[GFL]" );
        }
        
        if( _debug )
          printf("segmenter: %s (seqlen=%d (%d:%d),"
                 "nJumpsEst=%d, nJumps=%d, theta=%.2f)\n",
                 (gfl) ? "gflseg" : "heuristic",  shotlen, start_idx, end_idx,
                 njumps, (int)jump.size(), thrsh_gfl);
        
      }
      start_idx = end_idx = -1;
    }
  }
}

/*-----------------------------------------------------------------------*/
void VideoParser::sbd_heuristic(vector<double> v_diff, vector<int>& jump,
                                int njumps, int min_shot_len )
/*-----------------------------------------------------------------------*/
{
  vector<size_t> v_srt_idx;
  vector<double> v_srt_val; // contains sorted values
  hecate::sort( v_diff, v_srt_val, v_srt_idx );
  for(int i=(int)v_srt_val.size()-1; i>=0; i--) {
    bool add = true;
    if (((unsigned) v_srt_idx[i] + 1 < min_shot_len) ||
        ((unsigned) v_diff.size() - v_srt_idx[i] < min_shot_len)) {
      add = false;
    } else {
      for (size_t j = 0; j < jump.size(); j++) {
        int len = abs(jump[j] - (int) v_srt_idx[i]) + 1;
        if (len < min_shot_len) {
          add = false;
          break;
        }
      }
    }
    if( add ) {
      jump.push_back( (int)v_srt_idx[i] );
    }
    if( (int)jump.size() == njumps ) {
      break;
    }
  }
};

/*-----------------------------------------------------------------------*/
void VideoParser::extract_histo_features(int pyr_level, bool omit_filtered,
                                         int nbin_color, int nbin_edge_ori,
                                         int nbin_edge_mag)
/*-----------------------------------------------------------------------*/
{
  int npatches = 0;
  for(int l=0; l<pyr_level; l++)
    npatches += pow(4,l);
  
  int nbin_edge = nbin_edge_ori + nbin_edge_mag;
  Mat X_color_hist = Mat( npatches*3*nbin_color, _nfrm_given, CV_32F, Scalar(0,0,0)[0] );
  Mat X_edge_hist  = Mat( npatches*nbin_edge, _nfrm_given, CV_32F, Scalar(0,0,0)[0] );
  
#pragma omp parallel for
  for(int i=0; i<_nfrm_given; i++)
  {
    if( omit_filtered && !_v_frm_valid[i] ) continue;
    
    Mat color_hist;
    hecate::calc_pyr_color_hist( _v_frm_rgb[i], color_hist, nbin_color, pyr_level );
    color_hist.copyTo( X_color_hist.col(i) );
    
    Mat edge_hist;
    hecate::calc_pyr_edge_hist( _v_frm_gray[i], edge_hist, nbin_edge_ori, nbin_edge_mag, pyr_level );
    edge_hist.copyTo( X_edge_hist.col(i) );
    
  }
  
  // Transpose X_gray/edge_hist, X_hist = horzcat( X_color_hist, X_edge_hist )
  X_color_hist = X_color_hist.t();
  X_edge_hist = X_edge_hist.t();
  hconcat( X_color_hist, X_edge_hist, _X_feat );
}


/*-----------------------------------------------------------------------*/
void VideoParser::filter_redundant_and_obtain_subshots()
/*-----------------------------------------------------------------------*/
{
  if( _v_shot_ranges.empty() )
    update_shot_ranges();
  
  // Generate data matrix for kmeans, consider only valid frames
  int nfrm_valid = get_nfrm_valid();
  if( nfrm_valid==0 ) {
    return;
  }

  Mat km_data( nfrm_valid, _X_feat.cols, _X_feat.type() );
  vector<int> v_idxmap( nfrm_valid, 0 );
  
  int row=0;
  for(int i=0; i<_nfrm_given; i++)
  {
    if( _v_frm_valid[i] ) {
      _X_feat.row(i).copyTo(km_data.row(row));
      v_idxmap[row] = i;
      row++;
    }
  }
  
  //
  // Perform k-means
  int ncluster = min(nfrm_valid/2, (int)_v_shot_ranges.size());
  Mat km_lbl; // integer row vector; stores cluster IDs for every sample.
  Mat km_ctr; // One row per each cluster center.
  hecate::perform_kmeans( km_data, km_lbl, km_ctr, ncluster );
  
  //
  // convert km_lbl to v_frm_clusterid
  vector<int> v_frm_clusterid(_nfrm_given,-1);
  for(int i=0; i<km_lbl.rows; i++)
    v_frm_clusterid[ v_idxmap[i] ] = km_lbl.at<int>(i);
  
  //
  // Pick the most "still" frame from every sub-shot within each shot
  for(size_t shotid=0; shotid<_v_shot_ranges.size(); shotid++)
  {
    int sb0 = _v_shot_ranges[shotid].start;
    int sb1 = _v_shot_ranges[shotid].end;
    
    // Identify sub-shots (per kmeans result)
    int ssb0=-1, ssb1=-1, lbl=-1;
    for( int j=sb0; j<=sb1; j++ )
    {
      if( _v_frm_valid[j] ) {
        if( ssb0<0 ) { // enter a new sub-shot area
          ssb0 = j;
          lbl = v_frm_clusterid[j];
        }
        ssb1 = j;
      }
      if( ssb0>=0 && (v_frm_clusterid[j]!=lbl || j==sb1) )  // exit the current sub-shot area
      {
        int diff_min_idx=-1;
        double diff_min_val = numeric_limits<double>::max();
        for( int k=ssb0; k<=ssb1; k++ )
        {
          double diff_k = _X_diff.at<double>(k);
          if( diff_k < diff_min_val )
          {
            diff_min_idx = k;
            diff_min_val = diff_k;
          }
        }
        
        // Store subshot with keyframe index
        hecate::Range r(ssb0, ssb1);
        r.v_idx.push_back( diff_min_idx );
        _v_shot_ranges[shotid].v_idx.push_back( diff_min_idx );
        _v_shot_ranges[shotid].v_range.push_back( r );
        
        // Filter out redundant frames
        for( int k=ssb0; k<=ssb1; k++ )
          if( k!=diff_min_idx )
            mark_invalid(_v_frm_valid, _v_frm_log, k, "[Redundant]");
        
        // reset sub-shot area
        ssb0 = ssb1 = lbl = -1;
      }
    }
    
    // reset shot area
    sb0 = sb1 = -1;
  }
}


/*-----------------------------------------------------------------------*/
void VideoParser::update_shot_ranges( int min_shot_len )
/*-----------------------------------------------------------------------*/
{
  _v_shot_ranges.clear();
  
  int sb0=0, sb1=-1;
  for( int i=0; i<_nfrm_given; i++ )
  {
    if( _v_frm_valid[i] ) {
      // enter a new shot area
      if( sb0<0 ) sb0 = i;
      sb1 = i;
    }
    
    // exit the current shot area
    if( sb0>=0 && sb1>=0 && (!_v_frm_valid[i] || i+1==_nfrm_given) )
    {
      hecate::ShotRange r( sb0, sb1 );
      if( r.length()>min_shot_len ) {
        _v_shot_ranges.push_back( r );
      }
      else {
        for(int j=sb0; j<=sb1; j++ )
          mark_invalid( _v_frm_valid, _v_frm_log, j, "[Short]" );
      }
      sb0 = sb1 = -1;
    }
  }
}


/*-----------------------------------------------------------------------*/
int VideoParser::get_nfrm_valid()
/*-----------------------------------------------------------------------*/
{
  return (int) accumulate( _v_frm_valid.begin(), _v_frm_valid.end(), 0 );
}


/*-----------------------------------------------------------------------*/
void VideoParser::mark_invalid( vector<bool>& vec, int idx, int wnd_sz )
/*-----------------------------------------------------------------------*/
{
  int vec_len = (int)vec.size();
  for(int i=max(0,idx-wnd_sz); i<=min(vec_len-1,idx+wnd_sz); i++) {
    vec[i] = false;
  }
}

/*-----------------------------------------------------------------------*/
void VideoParser::mark_invalid( vector<bool>& vec, vector<string>& vec2,
                               int idx, const string msg, int wnd_sz )
/*-----------------------------------------------------------------------*/
{
  int vec_len = (int)vec.size();
  for(int i=max(0,idx-wnd_sz); i<=min(vec_len-1,idx+wnd_sz); i++) {
    vec[i] = false;
    vec2[i] = msg;
  }
}



/*-----------------------------------------------------------------------*/
void VideoParser::play_video_filtered( const string& in_video,
                                      int step_sz, int max_frm_len )
/*-----------------------------------------------------------------------*/
{
  printf("\nDebug mode: Displaying shot segmentation results...\n");
  
  // Debug visualization
  const char* WND = "DEBUG";
  namedWindow( WND ,1 );
  
  Mat frm;
  double debug_rsz_ratio = (double)max_frm_len/_video_width;
  
  VideoCapture vr( in_video );
  vr >> frm;
  resize( frm, frm, Size(), debug_rsz_ratio, debug_rsz_ratio, CV_INTER_LINEAR );
  vr.set( CV_CAP_PROP_POS_FRAMES, 0 );
  
  Size sz = frm.size();
  Mat frm_lr( sz.height, 2*sz.width, CV_8UC3 );
  Mat frm_l(frm_lr, Rect(0,0,sz.width,sz.height));
  Mat frm_r(frm_lr, Rect(sz.width,0,sz.width,sz.height));
  frm_lr.setTo(Scalar(0));
  
  for( int i=0; i<_nfrm_given; i++)
  {
    for( int j=0; j<step_sz; j++ ) {
      vr >> frm;
      if( frm.empty() ) break;
    }
    if( frm.empty() ) break;
    resize( frm, frm, Size(), debug_rsz_ratio, debug_rsz_ratio, CV_INTER_LINEAR );
    
    std::stringstream s;
    s << _v_frm_log[i];
    putText(frm, s.str(), Point2f(5,30), FONT_HERSHEY_PLAIN, 1.3, Scalar(255,0,255,255), 2);
    
    if( _v_frm_log[i].length()<2 ) {
      frm.copyTo( frm_l );
    }
    else {
      frm.copyTo( frm_r );
    }
    
    imshow( WND, frm_lr );
    if( waitKey(_video_fps)>=0 ) break;
  }
  vr.release();
  destroyWindow( WND );
}

