/*
 * Various helper functions for filesystem
 *
 * Copyright 2016 Yahoo Inc.
 * Licensed under the terms of the Apache 2.0 License.
 * See LICENSE file in the project root for terms.
 *
 * Developer: Yale Song (yalesong@yahoo-inc.com)
 */

#ifndef HECATE_FILE_HPP
#define HECATE_FILE_HPP

#include <stdio.h>
#include <string.h>

#include <iostream>
#include <fstream>
#include <algorithm>
#include <functional>
#include <cctype>
#include <locale>
#include <string>
#include <vector>

namespace hecate {
  
  // GET FILENAME WITHOUT EXTENSION AND TRAILING DIRECTORIES
  struct FileParts
  {
    std::string path;
    std::string name;
    std::string ext;
  };
  
  static inline FileParts fileparts(std::string filename)
  {
    int idx0 = filename.rfind("/");
    int idx1 = filename.rfind(".");
    
    if( idx1 == (int) std::string::npos )
      idx1 = filename.length();
    
    FileParts fp;
    fp.path = filename.substr(0,idx0+1);
    fp.name = filename.substr(idx0+1,idx1-idx0-1);
    fp.ext  = filename.substr(idx1);
    
    return fp;
  };
  
  static inline std::string get_dir( std::string filepath ) {
    hecate::FileParts fp = hecate::fileparts( filepath );
    std::string dir = fp.path;
    return dir;
  };
  
  static inline std::string escape_space( std::string s ) {
    std::string out;
    for( size_t i=0; i<s.size(); i++) {
      if( s[i] == ' ' )
        out += '\\';
      out += s[i];
    }
    return out;
  }
  
  static inline std::string get_filename( std::string filepath ) {
    hecate::FileParts fp = hecate::fileparts( filepath );
    std::string filename = fp.name;
    replace( filename.begin(), filename.end(), ' ', '_' );
    return filename;
  };
  
  // TRIM STRING. USEFUL FOR PROCESSING STRING FILENAMES
  // trim from start
  static inline std::string &ltrim(std::string &s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
    return s;
  };
  // trim from end
  static inline std::string &rtrim(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
    return s;
  };
  // trim from both ends
  static inline std::string &trim(std::string &s) {
    return ltrim(rtrim(s));
  };
  
  
  static inline char *trim(char *str)
  {
    size_t len = 0;
    char *frontp = str;
    char *endp = NULL;
    
    if( str == NULL ) { return NULL; }
    if( str[0] == '\0' ) { return str; }
    
    len = strlen(str);
    endp = str + len;
    
    /* Move the front and back pointers to address the first non-whitespace
     * characters from each end.
     */
    while( isspace(*frontp) ) { ++frontp; }
    if( endp != frontp )
    {
      while( isspace(*(--endp)) && endp != frontp ) {}
    }
    
    if( str + len - 1 != endp )
      *(endp + 1) = '\0';
    else if( frontp != str &&  endp == frontp )
      *str = '\0';
    
    /* Shift the string so that it starts at str so that if it's dynamically
     * allocated, we can still free it on the returned pointer.  Note the reuse
     * of endp to mean the front of the string buffer now.
     */
    endp = str;
    if( frontp != str )
    {
      while( *frontp ) { *endp++ = *frontp++; }
      *endp = '\0';
    }
    return str;
  };
  
  
  // CHECK IF FILE EXISTS
  static inline bool file_exists (const std::string& name) {
    if (FILE *file = fopen(name.c_str(), "r")) {
      fclose(file);
      return true;
    } else {
      return false;
    }
  };
  
  // READ TEXT FILE INTO A VECTOR
  static inline void read_textlist( const std::string& in,
                                   std::vector<std::string>& out )
  {
    std::string str;
    std::ifstream file( in );
    while( std::getline( file, str ) ) {
      str = hecate::trim(str);
      if( str.length()==0 ) continue;
      if( str.at(0)=='#' ) continue;
      out.push_back( str );
    }
  };
  
  // READ TEXT FILE CONTAINING A LIST OF FILEPATHS
  static inline void read_filelist( const std::string& in,
                                   std::vector<std::string>& out )
  {
    std::string str;
    std::ifstream file( in );
    while( std::getline( file, str ) ) {
      str = hecate::trim(str);
      if( str.length()==0 ) continue;
      if( str.at(0)=='#' ) continue;
      
      if( hecate::file_exists( str ) )
        out.push_back( str );
      else
        fprintf( stderr, "File doesn't exist: %s\n", str.c_str() );
      
    }
  };
  
  static inline void split_string( char* in,
                                  std::vector<std::string>& out,
                                  const char* delimiter)
  {
    out.clear();
    char *token = in;
    while( (token=strsep(&in,delimiter)) != NULL )
      out.push_back( hecate::trim(token) );
  };
  
  static inline std::string exec(const char* cmd) {
    FILE* pipe = popen(cmd, "r");
    if (!pipe) return "ERROR";
    char buffer[1024];
    std::string result = "";
    while (!feof(pipe)) {
      if (fgets(buffer, 1024, pipe) != NULL)
        result += buffer;
    }
    pclose(pipe);
    return result;
  };
  
  static inline std::string which(const char* bin) {
    char buf[1024];
    sprintf( buf, "echo `which %s`", bin );
    std::string result = hecate::exec(buf);
    result.erase( result.find_last_not_of(" \n\r\t")+1);
    return result;
  }
}
#endif

