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

#include "hecate/gflseg.hpp"

using namespace std;
using namespace cv;
using namespace hecate;

void Segmenter::gflseg(const Mat& X, vector<int>& jumps, int k, double theta)
{
  double epsilon = 1e-9;
  
  // group fused LASSO solver
  vector<int> jmp1;
  vector<double> lambda;
  gflars( X, jmp1, lambda, k, epsilon );
  
  // DP optimizer
  vector<vector<int> > jmp2;
  Mat rse;
  int kbest = 0;
  dpseg( X, jmp1, jmp2, rse, kbest, theta, -1 );
  
  // return the optimal solution
  jumps.resize( (int)jmp2[kbest].size(), 0 );
  for( size_t i=0; i<jmp2[kbest].size(); i++)
    jumps[i] = jmp2[kbest][i];
  //jumps[jmp2[kbest].size()] = X.rows-1; // last index
};


void Segmenter::dpseg(const Mat& X, const vector<int>& cp, vector<vector<int> >& jumps,
                      Mat& rse, int& kbest, double theta, int kmax )
{
  int n = X.rows; // the length of the signal X
  int p = X.cols; // the dimension of the signal X
  
  // Prevent the likelihood from diverging too much
  if( kmax<0 )
    kmax = min((int)cp.size(), (int)floor(n/10));
  
  vector<int>::iterator it;
  vector<int> cp_srt = cp;
  std::sort( cp_srt.begin(), cp_srt.end() );
  
  // Compute boundaries of the smallest intervals considered
  // b = sort(union([0],union([n],option.candidatechangepoints)));
  vector<int> b(2+cp.size());
  vector<int> edge = {-1, n-1};
  it = set_union( edge.begin(), edge.end(), cp_srt.begin(), cp_srt.end(), b.begin() );
  b.resize( it-b.begin() );
  std::sort( b.begin(), b.end() );
  
  
  // k is the number of such intervals
  int k = (int)b.size()-1;
  
  
  // Compute the k-by-k matrix J such that J(i,j) for i<=j is the RSE when
  // intervals i-to-j are merged.
  //
  // According to Picard,
  // "J(i,j) is the cost of the path connecting i-to-j in k steps (k segments)"
  // http://pbil.univ-lyon1.fr/members/fpicard/franckpicard_fichiers/pdf/aCGH-stat.pdf
  Mat J( k,   k, X.type(), Scalar(0,0,0)[0] );
  Mat S( n+1, p, X.type(), Scalar(0,0,0)[0] ); // cumsum of the rows of X
  Mat v( n+1, 1, X.type(), Scalar(0,0,0)[0] ); // cumsum of squared norm of the rows of X
  
  // S = [zeros(1,size(X,2));cumsum(X)]; % cumsum of the rows
  Mat Ssub = S( Rect(0,1,p,n) ); cumsum( X, Ssub );
  
  // v = [0;cumsum(sum(X.*X,2))]; % cumsum of squared norm of the rows
  Mat vsub = v( Rect(0,1,1,n) );
  cv::reduce( X.mul(X), vsub, 1, CV_REDUCE_SUM );
  cumsum( vsub, vsub );
  
  for( int i=0; i<k; i++ ) {
    for( int j=i; j<k; j++ ) {
      int Istart = b[i]+1;
      int Iend = b[j+1];
      // Regularization term: 1/length * ||s[end]-s[start]||^2_2
      double penalty = pow(cv::norm(S.row(Iend+1)-S.row(Istart)),2.0)/(Iend-Istart+1);
      J.at<double>(i,j) = v.at<double>(Iend+1)-v.at<double>(Istart)-penalty;
    }
  }
  
  //
  // DP recursion
  //
  
  // V(i,j) is the best RSE for segmenting intervals 1 to j with at most i-1 change-points
  Mat V( kmax+1, k, X.type(), Scalar(0,0,0)[0] );
  Mat jumpmat( kmax, k, CV_32S, Scalar(-1,0,0)[0] );
  
  // With no change-points, V(1,j) is just the precomputed RSE for intervals 1 to j
  J.row(0).copyTo( V.row(0) );
  
  // The recursive formula
  // NOTE: minidx/maxidx in minMaxIdx are very, very confusing to use!
  int minidx[2]; double minval;
  for( int ki=0; ki<kmax; ki++ ) {
    for( int j=ki+1; j<k; j++ ) {
      Mat tmp = V(Rect(ki,ki,j-ki,1)) + J(Rect(j,ki+1,1,j-ki)).t();
      minMaxIdx( tmp, &minval, 0, minidx, 0 ); // row matrix
      V.at<double>(ki+1,j) = minval;
      jumpmat.at<int>(ki,j) = minidx[1]+ki;
    }
  }
  
  
  // Optimal segmentations
  for( int ki=0; ki<kmax; ki++ )
  {
    vector<int> jump_ki(ki+1, 0);
    jump_ki[ki] = jumpmat.at<int>( ki,k-1 );
    for( int i=ki-1; i>=0; i-- )
      jump_ki[i] = jumpmat.at<int>(i,jump_ki[i+1]);
    jumps.push_back( jump_ki );
  }
  
  
  // Convert interval index back to the last position before the jump
  for( int ki=0; ki<kmax; ki++ )
    for( size_t i=0; i<jumps[ki].size(); i++ )
      jumps[ki][i] = b[ 1+jumps[ki][i] ];
  
  
  // RSE as a fuction of the number of change-points
  V.col(k-1).copyTo( rse );
  
  
  //
  // Based on DP table, find the optimal number of change-points
  //
  
  // log-likelihood is -n/2*(J-log(n)+1+log(2*pi));
  Mat JJ; cv::log( rse, JJ );
  
  // normalize
  // >> JJtild = (JJ(Km)-JJ) / (JJ(Km)-JJ(1)) * (Km-1) + 1;
  int Km = JJ.rows;
  Mat JJtild = (JJ.at<double>(Km-1)-JJ) / (JJ.at<double>(Km-1)-JJ.at<double>(0)) * (Km-1) + 1;
  
  
  // find the inflexion point
  // >> res.kbest = max(find(diff(diff(Jtild))>option.threshold))+1;
  Mat dJJtild(Km-1,1,X.type(),Scalar(0,0,0)[0]);
  Mat ddJJtild(Km-2,1,X.type(),Scalar(0,0,0)[0]);
  for(int i=0; i<Km-1; i++)
    dJJtild.at<double>(i) = JJtild.at<double>(i+1) - JJtild.at<double>(i);
  for(int i=0; i<Km-2; i++)
    ddJJtild.at<double>(i) = dJJtild.at<double>(i+1) - dJJtild.at<double>(i);
  
  kbest = 0;
  for( int i=0; i<Km-2; i++ )
    if( ddJJtild.at<double>(i)>theta && i>kbest )
      kbest = i;
};


void Segmenter::gflars(const Mat& X, vector<int>& jumps,
                       vector<double>& lambda, const int k, const double epsilon)
{
  int n = X.rows; // the length of the signal X
  int p = X.cols; // the dimension of the signal X
  
  jumps.clear();
  lambda.clear();
  
  // Default weight w(i) = sqrt(n/(i*(n-i)))
  Mat weights( n-1, 1, X.type() );
  for( int i=0; i<n-1; i++ )
    weights.at<double>(i)= sqrt((double)n/((i+1)*(n-i-1)));
  
  // Auxilary variable to use MATLAB-like sort function. The variable vjmps
  // should always be synced with vector<int> jumps.
  vector<int> vsrtval;    // contains sorted values
  vector<size_t> vsrtidx; // contains sorted indices
  
  //
  // Initialize cHat = W'*X
  Mat cHat;
  leftmultiplybyXt(X,weights,cHat);
  
  vector<int> A; // Active set indices (sorted in an ascending order)
  Mat cHatSub; // Used to access cHat(A,:)
  
  //
  // Main loop to find the successive jumps
  for( int iter=0; iter<k; iter++ )
  {
    // Compute row-wise norms of cHat
    // >> cHatSquareNorm = sum(cHat.^2,2);
    Mat cHatSquareNorm;
    norm2sq( cHat, cHatSquareNorm, 1 ); // 0:col, 1:row-wise
    
    // >> [bigcHat,besti]=max(cHatSquareNorm);
    int besti[2];
    double bigcHat;
    minMaxIdx( cHatSquareNorm, 0, &bigcHat, 0, besti ); // col matrix
    
    //
    // In the first iteration, we add the most correlated feature to the
    // active set. For the other iterations, this is already done at the
    // end of the previous iteration
    if( iter==0 ) {
      jumps.push_back( besti[0] );
    }
    
    // Resize active set vector and cHatSub matrix
    A.resize( iter+1, 0 );
    cHatSub = Mat( iter+1, p, X.type(), Scalar(0,0,0)[0] );
    
    
    //
    // Compute the descent direction W = inv(X(:,A)'*X(:,A))*cHat(A,:)
    Mat W; // size of (iter+1)-by-p
    
    // >> [A,I]=sort(res.jump(1:iter));
    hecate::sort( jumps, A, vsrtidx );
    get_submatrix_row( cHat, cHatSub, A );
    
    // >> w = leftmultiplybyinvXAtXA(n,A,cHat(A,:),weights);
    leftmultiplybyinvXAtXA(cHatSub,A,weights,n,W);
    
    // >> B = multiplyXtXbysparse(n,A,w,weights);
    Mat B;
    multiplyXtXbysparse(W,A,weights,n,B);
    
    //
    // Compute the descent step
    //   For each i we find the largest possible step alpha by solving:
    //      norm(cHat(i,:)-alpha*B(i,:)) = norm(cHat(j,:)-alpha*B(j,:))
    //      where j is in the active set.
    //   We write it as a second order polynomial
    //      a1(i)*alpha^2 - 2* a2(i)*alpha + a3(i)
    
    Mat a1,a2,a3;
    Mat a1sub, a2sub, a3sub;
    Mat tmp1, tmp2, tmp3;
    vector<int> subset;
    
    // >> a1 = bigcHat - sum(B.^2,2);
    cv::reduce( B.mul(B), a1, 1, CV_REDUCE_SUM );
    a1 = bigcHat - a1;
    
    // >> a2 = bigcHat - sum(B.*cHat,2);
    cv::reduce( B.mul(cHat), a2, 1, CV_REDUCE_SUM );
    a2 = bigcHat - a2;
    
    // >> a3 = bigcHat - cHatSquareNorm;
    a3 = bigcHat - cHatSquareNorm;
    
    //
    // Now we solve it
    // >> gammaTemp = zeros(2*(n-1),1);
    Mat gammaTemp( 2*(n-1), 1, X.type(), Scalar(0,0,0)[0] );
    
    // First, those where we really have a second-order polynomial
    // >> subset = find(a1 > EPSILON);
    for( int i=0; i<a1.rows; i++ )
      if( a1.at<double>(i)>epsilon )
        subset.push_back(i);
    
    if( !subset.empty() )
    {
      get_subvector( a1, a1sub, subset );
      get_subvector( a2, a2sub, subset );
      get_subvector( a3, a3sub, subset );
      
      // tmp1 = sqrt( a2(subset).^2 - a1(subset).*a3(subset) )
      cv::sqrt( a2sub.mul(a2sub)-a1sub.mul(a3sub), tmp1);
      
      // >> gammaTemp(subset)
      // = (a2(subset) + sqrt( a2(subset).^2 - a1(subset).*a3(subset) )) ./ a1(subset);
      cv::divide( a2sub + tmp1, a1sub, tmp2);
      for( size_t i=0; i<subset.size(); i++ )
        gammaTemp.at<double>(subset[i]) = tmp2.at<double>(i);
      
      // >> gammaTemp(subset+n-1)
      // = (a2(subset) - sqrt( a2(subset).^2 - a1(subset).*a3(subset) )) ./ a1(subset);
      cv::divide( a2sub - tmp1, a1sub, tmp2);
      for( size_t i=0; i<subset.size(); i++ )
        gammaTemp.at<double>(subset[i]+n-1) = tmp2.at<double>(i);
      
      subset.clear();
    }
    
    //
    // then those where the quadratic term vanishes and we have a
    // first-order polynomial
    // >> subset = find((a1 <= EPSILON) & (a2 > EPSILON));
    for( int i=0; i<a1.rows; i++ )
      if( a1.at<double>(i)<=epsilon && a2.at<double>(i)>epsilon )
        subset.push_back(i);
    
    if( !subset.empty() )
    {
      get_subvector( a1, a1sub, subset );
      get_subvector( a2, a2sub, subset );
      get_subvector( a3, a3sub, subset );
      
      // >> gammaTemp(subset)     = a3(subset) ./ (2*a2(subset));
      // >> gammaTemp(subset+n-1) = a3(subset) ./ (2*a2(subset));
      cv::divide( a3sub, 2*a2sub, tmp2 );
      for( size_t i=0; i<subset.size(); i++ )
        gammaTemp.at<double>(subset[i]) = gammaTemp.at<double>(subset[i]+n-1) = tmp2.at<double>(i);
      
      subset.clear();
    }
    
    
    //
    // Finally the active set should not be taken into account, as well as
    // those for which the computation gives dummy solutions
    // >> maxg=max(gammaTemp)+1;
    double maxg; minMaxIdx(gammaTemp, 0, &maxg); maxg+=1.0;
    
    // >> subset = find((a1 <= EPSILON) & (a2 <= EPSILON));
    for( int i=0; i<a1.rows; i++ )
      if( a1.at<double>(i)<=epsilon && a2.at<double>(i)<=epsilon )
        subset.push_back(i);
    
    if( !subset.empty() )
    {
      // >> gammaTemp(subset) = maxg;
      // >> gammaTemp(n+subset) = maxg;
      for( size_t i=0; i<subset.size(); i++ )
        gammaTemp.at<double>(subset[i]) = gammaTemp.at<double>(subset[i]+n) = maxg;
      subset.clear();
    }
    
    // >> gammaTemp(A) = maxg;
    // >> gammaTemp(n+A-1) = maxg;
    for( size_t i=0; i<A.size(); i++ )
      gammaTemp.at<double>(A[i]) = gammaTemp.at<double>(A[i]+n-1) = maxg;
    
    // >> gammaTemp(gammaTemp<=0)=maxg;
    // >> gammaTemp(imag(gammaTemp)~=0) = maxg;
    for( int i=0; i<gammaTemp.rows; i++ )
      if( gammaTemp.at<double>(i)<=0 || gammaTemp.at<double>(i)!=gammaTemp.at<double>(i) )
        gammaTemp.at<double>(i) = maxg;
    
    //
    // Now we can take the minimum
    // >> [gamma,nexttoadd]=min(gammaTemp);
    double gamma;
    int nexttoadd[2];
    minMaxIdx( gammaTemp, &gamma, 0, nexttoadd, 0); // col matrix
    
    //
    // Update
    // >> res.value{iter} = zeros(iter,p);
    // >> res.value{iter}(I,:) = gamma*w;
    // >> if iter>1
    // >>     res.value{iter}(1:(iter-1),:) = res.value{iter}(1:(iter-1),:) + res.value{iter-1};
    // >> end
    
    // >> res.lambda(iter)=sqrt(bigcHat);
    lambda.push_back( sqrt(bigcHat) );
    
    // >> if iter<k
    // >>     res.jump(iter+1) = 1+mod(nexttoadd-1,n-1);
    // >>     cHat = cHat-gamma*a;
    // >> end
    if( iter+1<k ) {
      jumps.push_back( nexttoadd[0]%(n-1) );
      cHat = cHat - gamma*B;
    }
  }
};


// X is n-by-p, w is (n-1)-by-1, R is (n-1)-by-p
void Segmenter::leftmultiplybyXt(const Mat& X, const Mat& w, Mat& R)
{
  int n = X.rows;
  int p = X.cols;
  
  R = Mat( n-1, p, X.type() );
  
  // R = ([1:n-1]'*U(end,:)/n - U(1:end-1,:)) .* w(:,ones(1,p));
  Mat U; cumsum( X, U );
  
  // llt = [1:n-1]'
  Mat llt( n-1, 1, X.type() );
  for( int r=0; r<n-1; r++ )
    llt.at<double>(r) = (double)(r+1);
  
  Mat lt = llt*U.row(n-1)/n - U(Rect(0,0,p,n-1));
  for( int c=0; c<p; c++ )
    R.col(c) = lt.col(c).mul(w);
};



// X is a-by-p, w is (n-1)-by-1, R is a-by-p, ind is a-by-1, where 1<=a<=(n-1)
void Segmenter::leftmultiplybyinvXAtXA(const Mat& X, const vector<int>& ind,
                                       const Mat& w, const int n, Mat& R)
{
  int a = X.rows;
  int p = X.cols;
  R = Mat( a, p, X.type(), Scalar(0,0,0)[0] );
  
  if( a>0 )
  {
    // >> u = diff([0;ind;n])
    // Note: we convert C++ index system (zero-base) to MATLAB (one-base)
    Mat u( a+1, 1, X.type() );
    u.at<double>(0) = ind[0]+1; // [2 0 1] becomes [3 1 2]
    u.at<double>(a) = n - (ind[a-1]+1);
    for( int i=1; i<a; i++ )
      u.at<double>(i) = ind[i] - ind[i-1];
    
    // >> val = val ./ w(ind,ones(1,p))
    Mat val;  X.copyTo(val);
    Mat wsub; get_subvector(w, wsub, ind);
    for( int c=0; c<p; c++ )
      cv::divide( val.col(c), wsub, val.col(c) );
    
    // >> delta = diff( [zeros(1,p); val; zeros(1,p)] ) ./ u(:,ones(1,p))
    Mat delta( a+1, p, X.type(), Scalar(0,0,0)[0] );
    delta.row(0) = val.row(0)+0.0;
    delta.row(a) = -val.row(a-1);
    for( int r=1; r<a; r++ )
      delta.row(r) = val.row(r) - val.row(r-1);
    for( int c=0; c<p; c++ )
      cv::divide( delta.col(c), u, delta.col(c) );
    
    // >> R = - diff( delta )
    for(int r=0; r<a; r++)
      R.row(r) = delta.row(r) - delta.row(r+1);
    
    // >> R = R ./ w(ind,ones(1,p))
    for(int c=0; c<p; c++)
      cv::divide( R.col(c), wsub, R.col(c) );
  }
};



void Segmenter::multiplyXtXbysparse(const Mat& X, const vector<int>& ind,
                                    const Mat& w, const int n, Mat& R)
{
  int a = X.rows;
  int p = X.cols;
  
  R = Mat( n-1, p, X.type(), Scalar(0,0,0)[0] );
  if( a>0 )
  {
    Mat wsub; get_subvector(w, wsub, ind);
    Mat val; X.copyTo(val);
    
    // First multiply beta by the weights
    // >> val = val .* w(ind,ones(1,p));
    for(int c=0; c<p; c++)
      val.col(c) = val.col(c).mul(wsub);
    
    //  compute the matrix s of increments of r
    // >> S = zeros(n-1,p);
    Mat S( n-1, p, X.type(), Scalar(0,0,0)[0] );
    
    // >> S(ind,:) = val;
    for( size_t i=0; i<ind.size(); i++ )
      val.row(i).copyTo(S.row(ind[i]));
    
    // >> S = flipud(cumsum(flipud(S)));
    flip(S,S,0); cumsum(S,S); flip(S,S,0);
    
    // >> u = ind' * val; // don't forget to add ones
    Mat indv( 1, a, X.type() );
    for( size_t i=0; i<ind.size(); i++ )
      indv.at<double>(i) = (double)(ind[i]+1);
    Mat u = indv * val;
    
    // >> S = S - u(ones(n-1,1),:)/n;
    for( int r=0; r<n-1; r++ )
      S.row(r) = S.row(r) - u/n;
    
    // then make the cumsum
    // >> R = cumsum(S);
    cumsum( S, R );
    
    // then multiply the rows by the weights
    // >> R = R .* w(:,ones(1,p));
    for( int c=0; c<p; c++ )
      R.col(c) = R.col(c).mul(w);
  }
};


