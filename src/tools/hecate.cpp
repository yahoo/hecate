/*
 * HECATE Yahoo Video Processing Library - Binary
 *
 * Copyright 2016 Yahoo Inc.
 * Licensed under the terms of the Apache 2.0 License.
 * See LICENSE file in the project root for terms.
 *
 * Developer: Yale Song (yalesong@yahoo-inc.com)
 */

#include "hecate/hecate.hpp"

using namespace std;
using namespace cv;

int main( int argc, char** argv )
{
  hecate_copyright();
  if( argc<3 )
    hecate_usage();
  
  // Read input params
  hecate_params opt;
  hecate_parse_params( argc, argv, opt );
  
  // Run VIDSUM
  vector<int> v_thumb_idx;
  vector<hecate::Range> v_gif_range;
  vector<hecate::Range> v_mov_range;
  run_hecate( opt, v_thumb_idx, v_gif_range, v_mov_range );
  
  // Print debugging info
  if( opt.debug ) {
    if( opt.jpg ) {
      printf("hecate: thumbnail indices: [ ");
      for(size_t i=0; i<v_thumb_idx.size(); i++)
        printf("%d ", v_thumb_idx[i]);
      printf("]\n");
    }
    if( opt.gif ) {
      printf("hecate: gif shot ranges:\n");
      for(int i=0; i<v_gif_range.size(); i++)
        v_gif_range[i].print();
    }
    if( opt.mov ) {
      printf("hecate: mov shot ranges:\n");
      for(int i=0; i<v_mov_range.size(); i++)
        v_mov_range[i].print();
    }
  }
  
  // Produce results
  if( opt.jpg ) {
    generate_thumbnails( opt, v_thumb_idx );
  }
  
  if( opt.gif ) {
    bool mov = opt.mov;
    opt.mov = false;
    generate_highlight_clips( opt, v_gif_range );
    opt.mov = mov;
  }
  
  if( opt.mov ) {
    bool gif = opt.gif;
    opt.gif = false;
    generate_highlight_clips( opt, v_mov_range );
    opt.gif = gif;
  }
}
