/*
 * 0/1 Knapsack solver
 *
 * Copyright 2016 Yahoo Inc.
 * Licensed under the terms of the Apache 2.0 License.
 * See LICENSE file in the project root for terms.
 *
 * Developer: Yale Song (yalesong@yahoo-inc.com)
 */

#ifndef HECATE_KNAPSACK_HPP
#define HECATE_KNAPSACK_HPP

#include <stdio.h>

namespace hecate {
  
  template <typename T>
  inline T Tmax(T a, T b) {return (a>b) ? a : b; };
  
  template <typename T>
  inline void solve_01knapsack( const vector<T>& value, const vector<int>& weight,
                               const int budget, vector<bool>& solution)
  {
    int n = (int)value.size();
    solution.resize(n, false);
    
    // recursion
    vector<vector<T> > V( n+1, vector<T>(budget+1,0) );
    for( int i=1; i<=n; i++ )
      for( int w=1; w<=budget; w++ )
        V[i][w] = ( weight[i-1]>w )
        ? V[i-1][w] : Tmax( V[i-1][w], V[i-1][w-weight[i-1]] + value[i-1] );
    
    // backtrack
    int w = budget;
    for( int i=n; i>0; i-- ) {
      if( V[i][w]!=V[i-1][w] && V[i][w]==V[i-1][w-weight[i-1]]+value[i-1] ) {
        solution[i-1] = true;
        w -= weight[i-1];
      }
    }
  }
  
}
#endif


