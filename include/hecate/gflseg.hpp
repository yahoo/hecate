/*
 * GFLSEG Group Fused LASSO Change Point Detection solver
 *
 * Copyright 2016 Yahoo Inc.
 * Licensed under the terms of the Apache 2.0 License.
 * See LICENSE file in the project root for terms.
 *
 * Developer: Yale Song (yalesong@yahoo-inc.com)
 *
 * Implementation is based on:
 *   Kevin Bleakley and Jean-Philippe Vert.
 *     "The group fused lasso for multiple change-point detection."
 *   arXiv preprint arXiv:1106.4199 (2011).
 */


#ifndef HECATE_GFLSEG_HPP
#define HECATE_GFLSEG_HPP

#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>

#if CV_MAJOR_VERSION > 2
#include <opencv2/core/core_c.h>
#endif

#include "hecate/sort.hpp"

using namespace std;
using namespace cv;

namespace hecate {
  
  class Segmenter {
  public:
    /*
     Automatic multiple change-point detection
     
     This function is a wrapper for multi-dimensional signal segmentation into at
     least k change-points by the LARS (function gflars), followed by change-point
     selection with a dynamic programming optimization (function dpseg).
     
     INPUT
     X      : n-by-p matrix to be segmented (n: # obs, p: # dim).
     k      : the number of change points to find in the first segmentation step
     theta  : stopping criteria, see dpseg. (default=0.1)
     forcek : k is set as the hard maximum number of change points
     
     OUTPUT
     jumps  : vector of change-points
     val    : matrix of the value on each interval
     updown : the up/down statistics on each interval
     */
    void gflseg(const Mat& X, vector<int>& jumps, int k, double theta=0.1);
    
    
    
    /*
     Segmentation of a multi-dimensional signal with dynamic programming
     http://pbil.univ-lyon1.fr/members/fpicard/franckpicard_fichiers/pdf/aCGH-stat.pdf
     Improvement is made by penalizing with the optimal segment length,
     lambda=sqrt(lstar-l) where lstar = n/k (typically chosen as 2-second long)
     
     INPUT
     X     : n-by-p matrix to be segmented (n: # obs, p: # dim).
     cp    : candidate change points (default: [0:n-2])
     kmax  : maximum number of change-points to test (default: length(cp))
     theta : stopping criteria. Typically chonse to be in the interval (0 0.5].
     The smaller the threshold, the higher the tendency to keep more
     breakpoints. The criteria is based on the method found in
     'Picard et al (2005)' "A statistical approach for array CGH data
     analysis" (BMC Bioinformatics). (default=0.5)
     
     OUTPUT
     jump  : j-by-1 vector of change-point positions for the i-th lambda value
     (j depends on lambda). i varies between 1 and kmax
     rse   : (kmax+1)-by-1 vector of residual squared error
     kbest : the number of selected change-points
     */
    void dpseg(const Mat& X, const vector<int>& cp, vector<vector<int> >& jumps,
               Mat& rse, int& kbest, double theta, int kmax=-1 );
    
    
    private:
    /*
     Segmentation of a multi-dimensional signal with the group fused LARS.
     
     INPUT
     X       : n-by-p matrix to be segmented (n: # obs, p: # dim).
     k       : the number of change points to find.
     epsilon : values smaller than epsilon are considered null.
     weights : (n-1)*1 vector of weights for the weighted graph fused LASSO.
     
     OUTPUT
     lambda     : estimated lambda values for each change-point
     jump       : successive change-point positions (1 x k)
     meansignal : mean signal per column (1 x p vector)
     */
    void gflars(const Mat& X, vector<int>& jumps, vector<double>& lambda,
                const int k, const double epsilon );
    
    
    /*
     Fast computation of Y' * X
     
     Compute R = Y'*X, where Y is the n-by-(n-1) design matrix for the weighted
     group fused lasso, with weights defined by the vector w, and X is a n-by-p
     data matrix. The computation is done in O(np).
     
     INPUT
     X : n-by-p matrix
     w : (n-1)-by-1 vector of weights
     
     OUTPUT
     R : (n-1)-by-p matrix equal to W' * X
     */
    void leftmultiplybyXt(const Mat& X, const Mat& w, Mat& R);
    
    
    /*
     Fast computation of inv(W'*W)*X
     
     Compute R = inv(W(:,ind)'*W(:,ind))*X, where W is the n-by-(n-1) design matrix
     
     INPUT
     X   : a-by-p matrix
     ind : a-by-1 vector of indices between 1 and n-1, sorted in an ascending order
     w   : (n-1)-by-1 vector of weights
     n   : the size of X is n-by-(n-1)
     
     OUTPUT
     R   : a-by-p matrix equal to inv(W'*W)*X
     */
    void leftmultiplybyinvXAtXA(const Mat& X, const vector<int>& ind,
                                const Mat& w, const int n, Mat& R);
    
    
    /*
     Fast computation of W' * W * X when X is sparse
     
     Compute R = W'*W*X, where X is a wor-sparse (n-1)-by-p matrix and W is the
     n-by-(n-1) design matrix
     
     INPUT
     X   : a-by-p matrix whose rows are the non-zero rows of the original X
     (same order as ind)
     ind : a-by-1 vector of indices of the non-zero rows of X (each in [1,n-1])
     w   : (n-1)-by-1 vector of weights
     n   : size of the problem
     
     OUTPUT
     R   : (n-1)-by-p matrix equal to W'*W*X
     */
    void multiplyXtXbysparse(const Mat& X, const vector<int>& ind,
                             const Mat& w, const int n, Mat& R);
    
    
    /* UTILITY FUNCTIONS */
    
    // dst = src(I)
    inline void get_subvector(const Mat& src, Mat& dst, const vector<int>& I)
    {
      dst = Mat( I.size(), 1, src.type() );
      for( size_t i=0; i<I.size(); i++ )
        dst.at<double>(i) = src.at<double>(I[i]);
    };
    
    // dst = src(I,:)
    inline void get_submatrix_row(const Mat& src, Mat& dst, const vector<int>& I)
    {
      dst = Mat( I.size(), src.cols, src.type() );
      for( size_t i=0; i<I.size(); i++ )
        src.row(I[i]).copyTo( dst.row(i) );
    };
    
    // dst = cumsum(src)
    inline void cumsum(const Mat& src, Mat& dst)
    {
      if( !dst.data )
        dst = Mat( src.rows, src.cols, src.type() );
      
      src.row(0).copyTo( dst.row(0) );
      for( int r=1; r<src.rows; r++)
        dst.row(r) = dst.row(r-1) + src.row(r);
    };
    
    // dst = sum(src.^2,2)
    inline void norm2sq( const Mat& src, Mat& dst, int dim)
    {
      if( dim==0 ) {// column-wise
        dst = Mat( 1, src.cols, src.type() );
        for(int c=0; c<src.cols; c++ )
          dst.at<double>(c) = cv::sum(src.col(c).mul(src.col(c)))[0];
      }
      else {
        dst = Mat( src.rows, 1, src.type() );
        for(int r=0; r<src.rows; r++)
          dst.at<double>(r) = cv::sum(src.row(r).mul(src.row(r)))[0];
      }
    }
    
    inline void check_nan(vector<double>& v) {
      for( size_t i=0; i<v.size(); i++ )
        if( v[i]!=v[i] )
          v[i]=0.0;
    };
    
    inline void print(Mat& X, string name)
    {
      cout << name << "_c=[" << endl;
      
      for(int r=0; r<X.rows; r++) {
        for( int c=0; c<X.cols; c++) {
          if( X.type()==CV_32F || X.type()==CV_64F )
            printf("%f ", X.at<double>(r,c));
          else
            printf("%d ", X.at<int>(r,c));
        }
        printf("\n");
      }
      printf("];\n\n");
    }
    
    inline void print(vector<int> v, string name)
    {
      cout << name << "_c=[";
      for(size_t i=0; i<v.size(); i++)
        printf("%d ", v[i]);
      printf("];\n\n");
    }
    
    inline void print(vector<double> v, string name)
    {
      cout << name << "_c=[";
      for(size_t i=0; i<v.size(); i++)
        printf("%f,", v[i]);
      printf("];\n\n");
    }
  };
  
}

#endif


