/*
 * Compute histogram of [x] using OpenCV
 *
 * Copyright 2016 Yahoo Inc.
 * Licensed under the terms of the Apache 2.0 License.
 * See LICENSE file in the project root for terms.
 *
 * Developer: Yale Song (yalesong@yahoo-inc.com)
 */

#ifndef HECATE_HIST_OPENCV_HPP
#define HECATE_HIST_OPENCV_HPP

#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>

#if CV_MAJOR_VERSION > 2
#include <opencv2/imgproc/types_c.h>
#endif

using namespace cv;
using namespace std;

namespace hecate {
  
  
  inline void print( const Mat x, const string filename )
  {
    ofstream myfile;
    myfile.open( filename );
    myfile << x;
    myfile.close();
  }
  
  inline void calc_gray_hist( const Mat& img, Mat& hist, int nbins=128 )
  {
    // Gray histogram
    float _range_val_gry[] = {0,256};
    const float* _range_gry = {_range_val_gry};
    
    calcHist( &img, 1, 0, Mat(), hist, 1, &nbins, &_range_gry, true, false);
    normalize( hist, hist );
  }
  
  inline void calc_color_hist( const Mat& img, Mat& hist, int nbins=128 )
  {
    // Color histogram
    float _range_val_clr[] = {0,256};
    const float* _range_clr = {_range_val_clr};
    
    // CV_BGR2HSV, CV_BGR2HLS, CV_BGR2YCrCb
    Mat img_cvt;
    cvtColor( img, img_cvt, CV_BGR2HSV );
    
    vector<Mat> planes;
    split( img_cvt, planes );
    
    if( !hist.data )
      hist = Mat( 3*nbins, 1, CV_32F, Scalar(0,0,0)[0] );
    
    Mat hist0( hist, Rect(0,0*nbins,1,nbins) );
    Mat hist1( hist, Rect(0,1*nbins,1,nbins) );
    Mat hist2( hist, Rect(0,2*nbins,1,nbins) );
    
    calcHist( &planes[0], 1, 0, Mat(), hist0, 1, &nbins, &_range_clr, true, false);
    calcHist( &planes[1], 1, 0, Mat(), hist1, 1, &nbins, &_range_clr, true, false);
    calcHist( &planes[2], 1, 0, Mat(), hist2, 1, &nbins, &_range_clr, true, false);
    
    normalize( hist0, hist0 );
    normalize( hist1, hist1 );
    normalize( hist2, hist2 );
  }
  
  
  inline void orientation( const Mat& Gx, const Mat& Gy, Mat& ori )
  {
    double alpha = 180.0 / 3.14159265358979323846264338327950288419716939937510;
    
    ori = Mat( Gx.rows, Gx.cols, CV_32F, Scalar(0,0,0)[0] );
    for(int r=0; r<Gx.rows; r++) {
      for(int c=0; c<Gx.cols; c++) {
        double deg = atan2(Gy.at<float>(r,c), Gx.at<float>(r,c)) * alpha;
        deg = (deg>=0) ? deg : deg + 360.0;
        deg = (deg<180.0) ? deg : deg - 180.0;
        ori.at<float>(r,c) = deg;
      }
    }
  }
  
  
  inline void calc_edge_hist( const Mat& Gx, const Mat& Gy, Mat& hist, int nbins_ori=16, int nbins_mag=16 )
  {
    // Edge histogram
    float _range_val_ori[] = {0,180};
    float _range_val_mag[] = {0,256};
    const float* _range_ori = {_range_val_ori};
    const float* _range_mag = {_range_val_mag};
    
    Mat ori, mag;
    orientation( Gx, Gy, ori );
    cv::magnitude( Gx, Gy, mag );
    
    if( !hist.data )
      hist = Mat( nbins_ori+nbins_mag, 1, CV_32F, Scalar(0,0,0)[0] );
    
    Mat hist_ori( hist, Rect(0,0,1,nbins_ori) );
    Mat hist_mag( hist, Rect(0,nbins_ori,1,nbins_mag) );
    
    calcHist( &ori, 1, 0, Mat(), hist_ori, 1, &nbins_ori, &_range_ori, true, false );
    calcHist( &mag, 1, 0, Mat(), hist_mag, 1, &nbins_mag, &_range_mag, true, false );
    
    normalize( hist_ori, hist_ori );
    normalize( hist_mag, hist_mag );
  }
  
  inline void calc_edge_hist( const Mat& img_gray, Mat& hist, int nbins_ori=16, int nbins_mag=16 )
  {
    Mat Gx, Gy;
    Scharr( img_gray, Gx, CV_32F, 1, 0 ); // ddepth, dx, dy
    Scharr( img_gray, Gy, CV_32F, 0, 1 );
    
    calc_edge_hist( Gx, Gy, hist, nbins_ori, nbins_mag );
  }
  
  
  inline void calc_pyr_gray_hist( const Mat& img, Mat& hist, int nbins=128, int level=2)
  {
    int w = img.cols;
    int h = img.rows;
    
    int npatches = 0;
    for(int i=0; i<level; i++)
      npatches += pow(4,i);
    
    int hist_sz = nbins;
    
    if( !hist.data )
      hist = Mat( hist_sz*npatches, 1, CV_32F, Scalar(0,0,0)[0] );
    
    int patch = 0;
    for(int l=0; l<level; l++) {
      for(int x=0; x<pow(2,l); x++) {
        for(int y=0; y<pow(2,l); y++) {
          int p_width  = floor( w/pow(2,l) );
          int p_height = floor( h/pow(2,l) );
          int p_x = x*p_width;
          int p_y = y*p_height;
          Mat patch_img( img, Rect(p_x,p_y,p_width,p_height) );
          Mat patch_hist( hist, Rect(0,patch*hist_sz,1,hist_sz) );
          calc_gray_hist( patch_img, patch_hist, nbins );
          patch++;
        }
      }
    }
  }
  
  
  inline void calc_pyr_color_hist( const Mat& img, Mat& hist, int nbins=128, int level=2)
  {
    int w = img.cols;
    int h = img.rows;
    
    int npatches = 0;
    for(int i=0; i<level; i++)
      npatches += pow(4,i);
    
    int hist_sz = 3 * nbins;
    
    if( !hist.data )
      hist = Mat( hist_sz*npatches, 1, CV_32F, Scalar(0,0,0)[0] );
    
    int patch = 0;
    for(int l=0; l<level; l++) {
      for(int x=0; x<pow(2,l); x++) {
        for(int y=0; y<pow(2,l); y++) {
          int p_width  = floor( w/pow(2,l) );
          int p_height = floor( h/pow(2,l) );
          int p_x = x*p_width;
          int p_y = y*p_height;
          Mat patch_img( img, Rect(p_x,p_y,p_width,p_height) );
          Mat patch_hist( hist, Rect(0,patch*hist_sz,1,hist_sz) );
          calc_color_hist( patch_img, patch_hist, nbins );
          patch++;
        }
      }
    }
  }
  
  
  inline void calc_pyr_edge_hist( const Mat& img, Mat& hist, int nbins_ori=16, int nbins_mag=16, int level=2)
  {
    int w = img.cols;
    int h = img.rows;
    
    int npatches = 0;
    for(int i=0; i<level; i++)
      npatches += pow(4,i);
    
    int hist_sz = nbins_ori + nbins_mag;
    
    if( !hist.data )
      hist = Mat( hist_sz*npatches, 1, CV_32F, Scalar(0,0,0)[0] );
    
    int patch = 0;
    for(int l=0; l<level; l++) {
      for(int x=0; x<pow(2,l); x++) {
        for(int y=0; y<pow(2,l); y++) {
          int p_width  = floor( w/pow(2,l) );
          int p_height = floor( h/pow(2,l) );
          int p_x = x*p_width;
          int p_y = y*p_height;
          Mat patch_img( img, Rect(p_x,p_y,p_width,p_height) );
          Mat patch_hist( hist, Rect(0,patch*hist_sz,1,hist_sz) );
          calc_edge_hist( patch_img, patch_hist, nbins_ori, nbins_mag );
          patch++;
        }
      }
    }
  }
  
  
  inline void calc_pyr_edge_hist( const Mat& Gx, const Mat& Gy, Mat& hist, int nbins_ori=16, int nbins_mag=16, int level=2)
  {
    int w = Gx.cols;
    int h = Gx.rows;
    
    int npatches = 0;
    for(int i=0; i<level; i++)
      npatches += pow(4,i);
    
    int hist_sz = nbins_ori + nbins_mag;
    
    if( !hist.data )
      hist = Mat( hist_sz*npatches, 1, CV_32F, Scalar(0,0,0)[0] );
    
    int patch = 0;
    for(int l=0; l<level; l++) {
      for(int x=0; x<pow(2,l); x++) {
        for(int y=0; y<pow(2,l); y++) {
          int p_width  = floor( w/pow(2,l) );
          int p_height = floor( h/pow(2,l) );
          int p_x = x*p_width;
          int p_y = y*p_height;
          Mat patch_Gx( Gx, Rect(p_x,p_y,p_width,p_height) );
          Mat patch_Gy( Gy, Rect(p_x,p_y,p_width,p_height) );
          Mat patch_hist( hist, Rect(0,patch*hist_sz,1,hist_sz) );
          calc_edge_hist( patch_Gx, patch_Gy, patch_hist, nbins_ori, nbins_mag );
          patch++;
        }
      }
    }
  }
  
}

#endif
