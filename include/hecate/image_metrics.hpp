/*
 * Various image metrics
 *
 * Copyright 2016 Yahoo Inc.
 * Licensed under the terms of the Apache 2.0 License.
 * See LICENSE file in the project root for terms.
 *
 * Developer: Yale Song (yalesong@yahoo-inc.com)
 */

#ifndef HECATE_IMAGE_METRICS_HPP
#define HECATE_IMAGE_METRICS_HPP

#include <math.h>
#include <stdio.h>
#include <cmath>
#include <string>
#include <vector>

#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>

#include <hecate/hist_opencv.hpp>

using namespace cv;
using namespace std;

namespace hecate {
  
  
#define VALIDATE(x) ( std::isnan(x) ) ? 0 : x
  
  /*------------------------------------------------------------------------
   Contrast is measured as the standard deviation of the pixel intensities.
   RMS contrast does not depend on the angular frequency content or the
   spatial distribution of contrast in the image
   https://en.wikipedia.org/wiki/Contrast_(vision)#RMS_contrast
   ------------------------------------------------------------------------*/
  inline double calc_rms_contrast( const Mat& gray_img )
  /*-----------------------------------------------------------------------*/
  
  {
    Mat I;
    gray_img.convertTo( I, CV_32FC1 );
    Scalar mu = cv::mean( I );
    I = I - mu[0];
    return VALIDATE(cv::norm(I)/sqrt(I.rows*I.cols));
  }
  
  /*------------------------------------------------------------------------
   Sharpness is measured as the sum of magnitude in frequency domain
   ------------------------------------------------------------------------*/
  inline double calc_sharpness( const Mat& gray_img )
  /*-----------------------------------------------------------------------*/
  {
    Mat img;
    gray_img.convertTo( img, CV_32FC1 );
    img *= 1./255;
    
    Mat dx, dy;
    Sobel( img, dx, img.type(), 1, 0, 3 );
    Sobel( img, dy, img.type(), 0, 1, 3 );
    magnitude( dx, dy, dx );
    
    int npixels = gray_img.rows * gray_img.cols;
    return VALIDATE(cv::sum(dx)[0] / npixels);
  }
  
  
  /*------------------------------------------------------------------------
   Brightness is measured as the relative luminance in colorimetric spaces
   https://en.wikipedia.org/wiki/Relative_luminance
   ------------------------------------------------------------------------*/
  inline double calc_brightness( const Mat& img )
  /*-----------------------------------------------------------------------*/
  {
    vector<Mat> bgr;
    split( img, bgr );
    for(size_t j=0; j<bgr.size(); j++)
      bgr[j].convertTo(bgr[j],CV_32F);
    
    // compute perceived brightness
    Mat img_pb = (0.2126*bgr[2] + 0.7152*bgr[1] + 0.0722*bgr[0])/255.0;
    return VALIDATE(mean( img_pb )[0]);
  }
  
  /*------------------------------------------------------------------------
   Uniformness is measured as the ratio between the 5% percentile pixel
   intensity histograms and the rest.
   ------------------------------------------------------------------------*/
  inline double calc_uniformity( const Mat& gray_img, int nbins=128 )
  /*-----------------------------------------------------------------------*/
  {
    double val = 0.0;
    
    Mat ghist;
    hecate::calc_gray_hist( gray_img, ghist, nbins );
    if( cv::sum(ghist)[0] == 0 ) {
      val = 1.0;
    }
    else {
      Mat hist_sorted;
      cv::sort( ghist, hist_sorted, CV_SORT_EVERY_COLUMN|CV_SORT_DESCENDING );
      hist_sorted /= cv::sum( hist_sorted )[0];
      for(int j=0; j<(int)(0.05*nbins); j++)
        val += (double) hist_sorted.at<float>(j);
    }
    return VALIDATE(val);
  }
  
  /*------------------------------------------------------------------------
   Symmetry is measured as the difference between edge orientation & magnitude
   histograms of the left and right halves.
   ------------------------------------------------------------------------*/
  inline double calc_asymmetry( const Mat& gray_img )
  /*-----------------------------------------------------------------------*/
  {
    Mat img_l( gray_img, Rect(0,0,gray_img.cols/2,gray_img.rows) );
    Mat img_r( gray_img, Rect(gray_img.cols/2,0,gray_img.cols/2,gray_img.rows) );
    
    Mat hist_l, hist_r;
    hecate::calc_edge_hist( img_l, hist_l );
    hecate::calc_edge_hist( img_r, hist_r );
    
    return VALIDATE((double)cv::norm( hist_l-hist_r ));
  }
  
  /*------------------------------------------------------------------------
   Third Saliency is measured as the residual saliency on the resized image,
   then collects the average saliency for each sub-window. It's an indicator
   of the informativeness of each window, but also of the general distribution
   of the content in the image.
   ------------------------------------------------------------------------*/
  inline vector<double> calc_third_saliency( const Mat& gray_img )
  {
    const int w = 3; // divide 128-by-128 image into 3-by-3 patches
    
    // 0-1 normalize & resize to 128-by-128
    Mat img;
    gray_img.convertTo( img, CV_32FC1 ); img *= 1./255;
    cv::resize( img, img, cv::Size(128,128), 0, 0, cv::INTER_AREA );
    
    //real is img image, im is zero
    Mat planes[]  = {Mat_<double>(img),               Mat::zeros(img.size(),CV_32FC1)};
    Mat planes2[] = {Mat::zeros(img.size(),CV_32FC1), Mat::zeros(img.size(),CV_32FC1)};
    Mat planesr[] = {Mat::zeros(img.size(),CV_32FC1), Mat::zeros(img.size(),CV_32FC1)};
    
    Mat complexI, magI, phasem;
    merge( planes, 2, complexI ); // Merges for FT
    dft( complexI, complexI );    //Forward DFT
    split( complexI, planes );    // planes[0] = Re(DFT(I), planes[1] = Im(DFT(I))
    magnitude( planes[0], planes[1], magI );    //DFT magnitude (amplitude)
    phase( planes[0], planes[1], phasem );      //DFT phase
    
    //log spectral amplitude
    cv::log( magI, magI );
    
    Mat smooth, kernel;
    smooth = Mat::zeros( img.size(), CV_32FC1 );
    kernel = Mat::zeros( Size(10,10), CV_32FC1 );
    kernel += Scalar::all( 0.01 ); //kernel for average filtering ( 0.01=1/(10*10) )
    
    //smoothed spectrum with average filtering
    cv::filter2D( magI, smooth, CV_32FC1, kernel, Point(-1,-1), 0, BORDER_REPLICATE );
    
    Mat diff;
    cv::subtract( magI, smooth, diff ); //spectral residual (log domain)
    cv::exp( diff, diff );              //spectral residual (real domain)
    
    //recover real and im part of the DFT after the residual
    polarToCart( diff, phasem, planes2[0], planes2[1] );
    
    //invert the DFT (we are back in the pixel domain! We have just created a Saliency Component)
    Mat result;
    merge( planes2, 2, complexI );
    dft( complexI, result, DFT_INVERSE+DFT_SCALE );
    split( result, planesr );                       // planes[0] = Re(DFT(I), planes[1] = Im(DFT(I))
    magnitude( planesr[0], planesr[1], result );	// sqrt(Re(DFT(I))^2 + Im(DFT(I))^2)
    
    normalize( result, result, 0, 2, cv::NORM_MINMAX );
    
    //compute the window dimensions
    float cx = result.cols/w;
    float cy = result.rows/w;
    vector<double> features;
    
    //local mean for each window
    for( int a=0; a<w; ++a ) {
      for( int b=0; b<w; ++b ) {
        Mat temp = result(Rect(cx*a, cy*b, cx, cy));
        double val = mean( temp )[0];
        features.push_back( VALIDATE(val) );
      }
    }
    
    return features;
  }
  /*-----------------------------------------------------------------------*/
  
  /*------------------------------------------------------------------------
   Entropy is measured -sum(p.*log2(p)), which is a statistical measure of
   randomness that can be used to characterize the texture of the imput image.
   http://www.mathworks.com/help/images/ref/entropy.html?refresh=true
   ------------------------------------------------------------------------*/
  inline double calc_entropy( const Mat& gray_img, int nbins=256 )
  /*-----------------------------------------------------------------------*/
  {
    double val=0.0;
    
    Mat ghist;
    hecate::calc_gray_hist( gray_img, ghist, nbins );
    for(int i=0; i<nbins; i++) {
      if( ghist.at<float>(i)!=0 )
        val -= ghist.at<float>(i) * log2(ghist.at<float>(i));
    }
    return VALIDATE(val);
  }
  
  
  /*------------------------------------------------------------------------
   Contrast balance is measured as the l2 distance between an original
   image and contrast-equalized image.
   ------------------------------------------------------------------------*/
  inline double calc_contrast_balance( const Mat& img )
  /*-----------------------------------------------------------------------*/
  {
    double res = 0;
    
    vector<Mat> hsv_planes;
    split( img, hsv_planes );
    
    //for each channel
    for (int pl=0; pl < 3; ++pl)
    {
      //equalize the histogram
      Mat img_hist_equalized;
      equalizeHist( hsv_planes[pl], img_hist_equalized );
      
      //compute the diff between original and equalized histograms
      res += cv::norm( hsv_planes[pl], img_hist_equalized ) /
      (double)(img.rows*img.cols);
    }
    
    //the larger the difference, the larger the contrast error. so we negate
    return VALIDATE(-res);
  }
  
  
  /*------------------------------------------------------------------------
   Exposure balance is measured as an absolute value of the skewness
   of the luminance histogram.
   ------------------------------------------------------------------------*/
  inline double calc_exposure_balance( const Mat& img )
  /*-----------------------------------------------------------------------*/
  {
    vector<Mat> hsv_planes;
    split( img, hsv_planes );
    
    /// Establish the number of bins
    const int histSize = 256;
    
    /// Set the ranges (for B,G,R)
    float range[] = { 0, 256 };
    const float *histRange = { range };
    
    cv::Mat hist;
    double res = 0;
    for (int pl=0; pl < 3; ++pl)
    {
      /// Compute the histograms:
      cv::calcHist( &hsv_planes[pl], 1, 0, cv::Mat(), hist, 1, &histSize, &histRange,
                   /*uniform: */true, /*accumulate: */false );
      cv::Scalar sump = cv::sum( hist );
      
      double mean = (double) sump[0] / hist.rows;
      
      cv::Mat prod;
      cv::multiply( (hist-mean), (hist-mean), prod);
      cv::Scalar sums = cv::sum( prod );
      
      cv::multiply( prod, (hist-mean), prod);
      cv::Scalar sumskew = cv::sum( prod );
      
      double result[2];
      result[0] = sqrt( (double) sums[0] / hist.rows );
      result[1] = (double) sumskew[0] / (hist.rows * result[0]*result[0]*result[0]);
      
      //absolute value of the skewness gives an idea of how much
      //the histogram of colors is incorrect
      res += std::abs( result[1] );
    }
    
    // the larger the skew, the larger the exposure error
    return VALIDATE(-res / 3);
  }
  
  /*------------------------------------------------------------------------
   JPEG quality is measured using the no-reference quality estimation of
   Z. Wang, H. R. Sheikh, and A. C. Bovik. No-reference perceptual quality
   assessment of jpeg compressed images. In ICIP, 2002.
   ------------------------------------------------------------------------*/
  inline double calc_jpeg_quality( const Mat& gray_img )
  /*-----------------------------------------------------------------------*/
  {
    // Convert to grayscale image with 64bit doubles
    cv::Mat img;
    gray_img.convertTo( img, CV_32FC1 );
    
    const int m = img.rows;
    const int n = img.cols;
    
    if (m<16 || n<16)
    {
      return -2.0;
    }
    
    // feature extraction: horizontal features
    Mat dh, dh1, dh2;
    dh1 = img( Rect(1, 0, n-1, m) );
    dh2 = img( Rect(0, 0, n-1, m) );
    dh = dh1 - dh2;
    
    double sum = 0;
    int count = 0;
    double sumz = 0;
    for (int i=0; i < m; ++i) {
      for (int j=0; j < n-2; ++j) {
        if ((j+1)%8==0 && j>0 && (j+1)<8*floor(n/8))
        {
          sum += std::abs( dh.at<double>(i,j) );
          count++;
        }
        double signval = copysign( 1.0, dh.at<double>(i,j) ) *
        copysign( 1.0, dh.at<double>(i,j+1) );
        if (signval < 0)
          sumz += 1;
      }
    }
    
    double bh = sum / count;
    double ah = (8.0 * mean( cv::abs(dh) )[0] - bh) / 7;
    double zh = sumz / (m * (n-2));
    
    // feature extraction: vertical features
    Mat dv1, dv2, dv;
    dv1 = img( Rect(0, 1, n, m-1) );
    dv2 = img( Rect(0, 0, n, m-1) );
    dv = dv1 - dv2;
    
    sum = 0;
    count = 0;
    sumz = 0;
    for (int i=0; i < m-2; ++i) {
      for (int j=0; j < n; ++j) {
        if ((i+1)%8==0 && i>0 && (i+1)<8*floor(m/8))
        {
          sum += std::abs( dv.at<double>(i,j) );
          count++;
        }
        double signval = copysign( 1.0, dv.at<double>(i,j) ) *
        copysign( 1.0, dv.at<double>(i,j+1) );
        if (signval < 0)
          sumz += 1;
      }
    }
    
    double bv = sum / count;
    double av = (8.0 * mean( cv::abs(dv) )[0] - bv) / 7;
    double zv = sumz / (n * (m-2));
    
    //combined features
    double B = (bh + bv) / 2;
    double A = (ah + av) / 2;
    double Z = (zh + zv) / 2;
    
    // Quality Prediction
    const double alpha = -245.8909;
    const double beta = 261.9373;
    const double gamma1 = -239.8886;
    const double gamma2 = 160.1664;
    const double gamma3 = 64.2859;
    
    double score = alpha + beta *
    std::pow( B, (gamma1/10000) ) *
    std::pow( A, (gamma2/10000) ) *
    std::pow( Z, (gamma3/10000) );
    
    return VALIDATE(score);
  }
  
}
#endif
