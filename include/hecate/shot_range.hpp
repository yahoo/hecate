/*
 * Shot range definition
 *
 * Copyright 2016 Yahoo Inc.
 * Licensed under the terms of the Apache 2.0 License.
 * See LICENSE file in the project root for terms.
 *
 * Developer: Yale Song (yalesong@yahoo-inc.com)
 */

#ifndef HECATE_SHOT_RANGE_HPP
#define HECATE_SHOT_RANGE_HPP

#include <string>
#include <vector>

#include <hecate/sort.hpp>

namespace hecate {
  
  /*
   */
  class Range {
  public:
    int start;
    int end;
    vector<int> v_idx; // store indices (keyframe index)
    
  public:
    Range(int s, int e): start(s), end(e) {};
    
    inline int length() const { return end-start+1; };
    
    inline void print() const {
      printf("range(%d:%d) (%d) [", start, end, length() );
      for(size_t i=0; i<v_idx.size(); i++) {
        printf("%d", v_idx[i]);
        if( i+1<v_idx.size() )
          printf(",");
      }
      printf("]\n");
    };
  };
  
  /*
   */
  class ShotRange: public Range {
  public:
    vector<Range> v_range;
    
  public:
    ShotRange(int s, int e): Range(s,e) {};
    
    inline void print() const {
      printf("shot_range(%d:%d) (%d) [", start, end, length() );
      for(size_t i=0; i<v_idx.size(); i++) {
        printf("%d", v_idx[i]);
        if( i+1<v_idx.size() )
          printf(",");
      }
      printf("]\n");
      
      for(size_t i=0; i<v_range.size(); i++) {
        printf("  sub_");
        v_range[i].print();
      }
    };
  };
  
  /*
   */
  class Tag {
  public:
    string label;
    double score;
    
  public:
    Tag(string& l, double s): label(l), score(s) {};
    
    inline void print() const {
      printf("(%s,%f)", label.c_str(), score);
    }
  };
  
  /*
   */
  class RangeTag: public Range {
  public:
    vector<Tag> v_tag;
  
  public:
    RangeTag(int s, int e): Range(s,e) {};
    RangeTag(Range r): Range(r.start,r.end) {
      v_idx=r.v_idx;
    };
    
    inline void print() const {
      printf("range(%d:%d) (%d) [", start, end, length() );
      for(size_t i=0; i<v_idx.size(); i++) {
        printf("%d", v_idx[i]);
        if( i+1<v_idx.size() )
          printf(",");
      }
      printf("], tags: ");
      for(size_t i=0; i<v_tag.size(); i++)
        v_tag[i].print();
      printf("\n");
    };
    
    inline void sort() {
      vector<double> v_scores;
      for(size_t i=0; i<v_tag.size(); i++)
        v_scores.push_back( v_tag[i].score );
      vector<double> v_srtval;
      vector<size_t> v_srtidx;
      hecate::sort( v_scores, v_srtval, v_srtidx );
      vector<Tag> v_tag_new;
      for(size_t i=0; i<v_srtidx.size(); i++)
        v_tag_new.push_back( v_tag[v_srtidx[v_srtidx.size()-1-i]] );
      v_tag = v_tag_new;
    }
  };
  
  /*
   */
  class ShotRangeTag: public Range {
  public:
    vector<RangeTag> v_range_tag;
    
  public:
    ShotRangeTag(int s, int e): Range(s,e) {};
    ShotRangeTag(ShotRange& sr): Range(sr.start,sr.end) {
      v_idx = sr.v_idx;
      for(size_t i=0; i<sr.v_range.size(); i++) {
        RangeTag rt( sr.v_range[i] );
        v_range_tag.push_back( rt );
      }
    };
    
    inline void print() {
      printf("shot_range_tag(%d:%d) (%d) [", start, end, length() );
      for(size_t i=0; i<v_idx.size(); i++) {
        printf("%d", v_idx[i]);
        if( i+1<v_idx.size() )
          printf(",");
      }
      printf("]\n");
      
      for(size_t i=0; i<v_range_tag.size(); i++) {
        printf("  sub_");
        v_range_tag[i].print();
      }
    };
    
    inline void sort() {
      for(size_t i=0; i<v_range_tag.size(); i++)
        v_range_tag[i].sort();
    };
  };
  
}
#endif


