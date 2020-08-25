/*
 * Gap statistics
 *
 * Copyright 2016 Yahoo Inc.
 * Licensed under the terms of the Apache 2.0 License.
 * See LICENSE file in the project root for terms.
 *
 * Developer: Yale Song (yalesong@yahoo-inc.com)
 */


#ifndef HECATE_GAPSTAT_HPP
#define HECATE_GAPSTAT_HPP

#include <stdio.h>
#include <stdlib.h>

#include <iostream>   // for standard I/O
#include <limits>
#include <numeric>

#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>

#if CV_MAJOR_VERSION > 2
#include <opencv2/core/core_c.h>
#include <opencv2/imgproc/types_c.h>
#endif

using namespace std;
using namespace cv;

namespace hecate {
  
  // Compute mean-squared error (mse)
  inline double calc_mse( const Mat& km_data, const Mat& km_lbl, const Mat& km_ctr )
  {
    double compactness = 0;
    
    int nsamples = km_lbl.rows;
    int nclusters  = km_ctr.rows;
    
    for( int k=0; k<nclusters; k++ ) {
      int nk = 0;
      double err = 0;
      for( int i=0; i<nsamples; i++ ) {
        if( km_lbl.at<int>(i)==k ) {
          nk++;
          err += pow(cv::norm(km_data.row(i)-km_ctr.row(k)),2.0);
        }
      }
      compactness += 0.5 * err / nk;
    }
    
    return compactness;
  }
  
  inline void find_bounds( const Mat& km_data, Mat& bounds )
  {
    int ndim = km_data.cols;
    bounds = Mat(2,ndim,CV_64F);
    
    double minval, maxval;
    for( int d=0; d<ndim; d++ ){
      minMaxLoc( km_data.col(d), &minval, &maxval );
      bounds.at<double>(0,d) = minval;
      bounds.at<double>(1,d) = maxval;
    }
  }
  
  inline void randu_bound( Mat& X, const Mat& bounds, int nsamples, int ndims )
  {
    X = Mat(nsamples, ndims, CV_32F);
    randu( X, Scalar::all(0.0), Scalar::all(1.0) );
    for(int c=0; c<X.cols; c++) {
      double lb = bounds.at<double>(0,c);
      double ub = bounds.at<double>(1,c);
      X.col(c) = (ub-lb)*X.col(c) + lb;
    }
  }
  
  //fix error when there is only 1 data point: http://docs.opencv.org/2.4/modules/core/doc/clustering.html
  inline void perform_kmeans(const Mat& km_data, Mat& km_lbl, Mat& km_ctr, int ncluster,
                             int km_attempts=1, int km_max_cnt=1000, double km_eps=0.0001)
  {
      if(km_data.rows==1){
          km_lbl = Mat::zeros(1,1, km_lbl.type());
          cv::reduce( km_data, km_ctr, 0, CV_REDUCE_AVG );
      }
      else{
        int km_k = min(ncluster, km_data.rows);
        TermCriteria km_opt = TermCriteria(CV_TERMCRIT_ITER|CV_TERMCRIT_EPS, km_max_cnt, km_eps);
        kmeans( km_data, km_k, km_lbl, km_opt, km_attempts, KMEANS_PP_CENTERS, km_ctr );
      }
}
  
  inline int perform_kmeans_gs( const Mat& km_data, Mat& km_lbl, Mat& km_ctr,
                               vector<int> K, int B=10, int N=500 )
  {
    int ndims = km_data.cols;
    
    // Generate null reference dataset
    vector<Mat> v_Xb;
    Mat bounds; find_bounds( km_data, bounds );
    for(int b=0; b<B; b++) {
      Mat Xb;
      randu_bound( Xb, bounds, N, ndims );
      v_Xb.push_back( Xb );
    }
    
    // Compute Gap Statistics
    vector<double> gap(K.size(), 0.0);    // gap statistics = E[logWk] - logWk
    vector<double> ElogWk(K.size(), 0.0); // compactness of the null reference distribution
    vector<double> logWk(K.size(), 0.0);  // compactness of the given data
    vector<double> Sk(K.size(), 0.0);   // standard deviations
    
    for( size_t i=0; i<K.size(); i++ )
    {
      Mat lbl, ctr;
      perform_kmeans( km_data, lbl, ctr, K[i] );
      logWk[i] = log(calc_mse( km_data, lbl, ctr ));
      
      ElogWk[i] = 0.0;
      vector<double> logWkb(B, 0.0);
      for(int b=0; b<B; b++)
      {
        Mat lbl_b, ctr_b;
        perform_kmeans( v_Xb[b], lbl_b, ctr_b, K[i], 100, 1000, 0.001 );
        logWkb[b] = log(calc_mse( v_Xb[b], lbl_b, ctr_b ));
        ElogWk[i] += logWkb[b];
      }
      ElogWk[i] /= B;
      
      Sk[i] = 0.0;
      for(int b=0; b<B; b++)
        Sk[i] += (logWkb[b]-ElogWk[i])*(logWkb[b]-ElogWk[i]);
      Sk[i] = sqrt(1.0+1.0/B) * sqrt( Sk[i] / B );
      
      gap[i] = ElogWk[i] - logWk[i];
      printf("\tgapstat: k=%d, logWk=%f, ElogWk=%f, Sk=%f, Gap=%f\n",
             K[i], logWk[i], ElogWk[i], Sk[i], gap[i]);
    }
    
    // find the smallest k such that gap(k) >= gap(k+1) - S[k+1]
    size_t kstar=0;
    for(size_t i=0; i<K.size()-1; i++) {
      if( gap[i] >= gap[i+1] - Sk[i+1] ) {
        kstar = i;
        break;
      }
    }
    printf("\tgapstat: Optimal k=%d [%d:%d]\n", K[kstar], K[0], K[K.size()-1]);
    
    // return results
    perform_kmeans( km_data, km_lbl, km_ctr, K[kstar] );
    return K[kstar];
    
  }
}

#endif

