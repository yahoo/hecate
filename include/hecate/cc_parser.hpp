/*
 * Video caption parser (VTT and TTML)
 *
 * Copyright 2016 Yahoo Inc.
 * Licensed under the terms of the Apache 2.0 License.
 * See LICENSE file in the project root for terms.
 *
 * Developer: Yale Song (yalesong@yahoo-inc.com)
 */

#ifndef HECATE_CC_PARSER_HPP
#define HECATE_CC_PARSER_HPP

#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cstring>
#include <algorithm>

#include "hecate/time.hpp"

using namespace std;

namespace hecate {
  
  struct tcc {
    double start; // milli-seconds
    double end;   // milli-seconds
    string text;
  };
  
  inline double time_str2num( const string str_time )
  {
    int hh,mm,ss,ms;
    
    hh = stoi( str_time.substr(0,2) );
    mm = stoi( str_time.substr(3,2) );
    ss = stoi( str_time.substr(6,2) );
    ms = stoi( str_time.substr(9,3) );
    
    return 60*60*hh + 60*mm + ss + (double)ms/1000.0;
  }
  
  
  inline string truncate_caption( const string src, const int maxlen )
  {
    string dst = src;
    while( (int)dst.length()>maxlen )
    {
      size_t pos = dst.find_last_of(" ");
      dst = dst.substr(0,pos);
    }
    return dst;
  }
  
  inline string clean_caption( const string src )
  {
    string dst = src;
    size_t pos, pos2;
    
    // Remove special characters
    char chars[] = ".,;:@#$%^&*\"><_+=~`/";
    for(unsigned i=0; i<strlen(chars); i++)
      dst.erase( remove(dst.begin(), dst.end(), chars[i]), dst.end() );
    
    // Remove leading & trailing spaces
    dst.erase(dst.begin(), find_if(dst.begin(), dst.end(), bind1st(not_equal_to<char>(), ' ')));
    dst.erase(find_if(dst.rbegin(), dst.rend(), bind1st(not_equal_to<char>(), ' ')).base(), dst.end());
    
    // Remove leadning words
    vector<string> lwords = {"And ","But ","Or ","and ","Then ","then "};
    for(size_t i=0; i<lwords.size(); i++)
    {
      pos = dst.find(lwords[i]);
      if( pos==0 ) {
        //printf("<<[%s>> [%s] -> [%s]\n", lwords[i].c_str(), dst.c_str(), dst.substr(lwords[i].length(),string::npos).c_str());
        dst = dst.substr(lwords[i].length(),string::npos);
      }
    }
    
    // Remove trailing words
    dst = dst + " ";
    vector<string> words = {" and "," or "," but "," however "," to "," of "," so ",
      " because "," who "," when "," where "," what "," how ", " why ",
      " while "," which "," with "," the "};
    for(size_t i=0; i<words.size(); i++)
    {
      pos = dst.substr(0,dst.find_last_of(" ")).find_last_of(" ");
      if( pos==string::npos )
        break;
      pos2 = dst.substr(pos,string::npos).find(words[i]);
      if( pos2!=string::npos ) {
        //printf("<<%s]>> [%s] -> [%s]\n", words[i].c_str(), dst.c_str(), dst.substr(0,pos).c_str());
        dst = dst.substr(0,pos) + " ";
      }
    }
    
    // Remove annotation inside a bracket, e.g., [UNKNOWN], [LAUGH], [MUSIC], etc.
    pos = dst.find("[");
    while( pos!=string::npos )
    {
      pos2 = dst.find("]");
      dst = dst.substr(0,pos) + dst.substr(pos2+1,string::npos);
      pos = dst.find("[");
    }
    
    // Remove double ??
    pos = dst.find("??");
    while( pos!=string::npos )
    {
      dst = dst.substr(0,pos) + dst.substr(pos+2,string::npos);
      pos = dst.find("??");
    }
    
    // Remove double space
    pos = dst.find("  ");
    while( pos!=string::npos )
    {
      dst = dst.substr(0,pos) + dst.substr(pos+1,string::npos);
      pos = dst.find("  ");
    }
    
    // Remove leading & trailing spaces
    dst.erase(dst.begin(), find_if(dst.begin(), dst.end(), bind1st(not_equal_to<char>(), ' ')));
    dst.erase(find_if(dst.rbegin(), dst.rend(), bind1st(not_equal_to<char>(), ' ')).base(), dst.end());
    
    
    return dst;
  }
  
  
  
  inline void print_closed_caption( vector<tcc>& cc )
  {
    for(size_t i=0; i<cc.size(); i++ )
      printf("%.2f  -->  %.2f: [%s]\n",cc[i].start,cc[i].end,cc[i].text.c_str());
  }
  
  
  
  inline void parse_vtt( const string filename, vector<tcc>& vtt )
  {
    string line;
    
    ifstream fid;
    fid.open( filename );
    getline( fid, line );
    while( fid )
    {
      size_t pos = line.find(" --> ");
      if( pos!=string::npos )
      {
        tcc cc;
        
        // time
        cc.start = time_str2num( line.substr( 0, pos ) );
        cc.end = time_str2num( line.substr( pos+5, string::npos) );
        
        // text
        getline( fid, line );
        line = clean_caption( line );
        
        pos = line.find("???");
        if( pos==string::npos && line.length()>3 ) {
          cc.text = line;
          vtt.push_back( cc );
        }
      }
      getline( fid, line );
    }
  }
  
  inline void parse_ttml( const string filename, vector<tcc>& ttml )
  {
    size_t pos1, pos2, pos3;
    string line, item;
    
    string p_start = "<p ";
    string p_end   = "</p>";
    
    ifstream fid;
    fid.open( filename );
    getline( fid, line );
    
    // Check the language
    pos1 = line.find("<language>en-US</language>");
    if( pos1==string::npos )
      return;
    
    // Parse TTML
    while( fid )
    {
      pos1 = line.find(p_start);
      while( pos1!=string::npos )
      {
        pos2 = line.find(p_end);
        item = line.substr(pos1, pos2-pos1);
        line = line.substr( pos2+3, string::npos );
        
        pos1 = item.find("begin=");
        pos2 = item.find("end=");
        pos3 = item.find("\">");
        
        tcc cc;
        cc.start = time_str2num( item.substr(pos1+7,12) );
        cc.end   = time_str2num( item.substr(pos2+5,12) );
        cc.text  = clean_caption( item.substr(pos3+2,string::npos) );
        ttml.push_back( cc );
        
        pos1 = line.find(p_start);
      }
      getline( fid, line );
    }
    
    //print_closed_caption( ttml );
  }
  
  //
  inline string encode_vtt( int index, float start_sec, float end_sec,
                           vector<string>& v_msg )
  {
    char buf[512];
    string start_sec_str = hecate::second2string( start_sec, "hh:mm:ss.mss" );
    string end_sec_str = hecate::second2string( end_sec, "hh:mm:ss.mss" );
    sprintf( buf, "%d\n%s --> %s\n", index, start_sec_str.c_str(),
            end_sec_str.c_str());

    for(size_t i=0; i<v_msg.size(); i++) {
      sprintf( buf, "%s- %s\n", buf, v_msg[i].c_str());
    }
    
    string ret = buf;
    return ret;
  }
  
}
#endif

