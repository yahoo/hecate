# HECATE
Hecate [hek-uh-tee] is a video processing library that auto-magically generates thumbnails, animated GIFs, and video summaries from videos. This library is developed and maintained by Yahoo Research, New York.

The source code is Copyright 2016 Yahoo Inc. and is licensed under the terms of the Apache 2.0 License. See the [LICENSE](https://github.com/yahoo/hecate/blob/master/LICENSE) in the project root file for terms.

The technology behind this library is based on our research work. If you find this library useful in your work, we ask you to cite our research paper:
```
"To Click or Not To Click: Automatic Selection of Beautiful Thumbnails from Videos."
Yale Song, Miriam Redi, Jordi Vallmitjana, Alejandro Jaimes, 
Proceedings of the 25th ACM International on Conference on Information and Knowledge Management, CIKM 2016
```

## Installation
Hecate has one dependency: [OpenCV library](https://github.com/opencv/opencv) with an [FFMPEG](https://github.com/FFmpeg/FFmpeg) support. You will need to install the library properly before trying out Hecate!

Once you install the dependenct library correctly, follow the instruction below:
```
$ git clone https://github.com/yahoo/hecate.git
$ cd hecate
$ vim Makefile.config
 - Set INCLUDE_DIRS and LIBRARY_DIRS to where your 
   opencv library is installed. Usually under /usr/local.
 - If your OpenCV version is 2.4.x, comment out the line 
   OPENCV_VERSION := 3
 - Save and exit
$ make all
$ make distribute
```

Once you've successfully compiled hecate, it will generate a binary executable under `distribute/bin/`. Run the following command to check if everything works properly:
```
$ ./distribute/bin/hecate

-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
 HECATE Copyright 2016 Yahoo Inc.
   Licensed under the terms of the Apache 2.0 License.
   Developed by : Yale Song (yalesong@yahoo-inc.com)
   Built on  : 11:46:03 Aug 11 2016
-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
USAGE: hecate -i infile [options]

  -i  --in_video      (string)    Input video file
  -o  --out_dir       (string)    Output directory (./output)
  -s  --step          (int)       Frame subsampling step size (1)
  -n  --njpg          (int)       Number of thumbnails to be generated (5)
  -q  --ngif          (int)       Number of GIFs to be generated (5)
  -r  --lmov          (int)       Length of video summary to be generated (in seconds) (15)
  -u  --jpg_width_px  (int)       Pixel width of thumbnail images (360)
  -v  --gif_width_px  (int)       Pixel width of animated GIFs (360)
  -w  --mov_width_px  (int)       Pixel width of summary video (360)
  --generate_jpg                  Generate thumbnail images
  --generate_gif                  Generate animated GIFs
  --generate_mov                  Generate a summary video
  --generate_gifsum               Generate animated GIFs summary
  --generate_gifall               Generate all possible animated GIFs
  --print_shot_info               Print shot boundary detection results
  --print_keyfrm_info             Print keyframe indices
```

Congratulations! You have successfully installed hecate!


## Get started
In order to get started, we will need a video file to play with. In this example, we will use the video ["The Spirit of '43" by Walt Disney](https://archive.org/details/TheSpiritOf43_56) from [The Internet Archive](https://archive.org). 

Let's download the video and save it as `examples/video.mp4`:
```
$ wget https://archive.org/download/TheSpiritOf43_56/The_Spirit_of__43_512kb.mp4 \
  --output-document examples/video.mp4 --no-check-certificate
```

Hecate provides three main functionalities through a binary executable `hecate`: Thumbnail extraction, GIF generation, and video summarization. There are various other functionalities the library provides, such as shot boundary detection and keyframe extraction. 

We will explain each case below.

### Shot boundary detection and keyframe extraction
Shot boundary detection and keyframe extraction are often the first steps towards various video processing methods. With Hecate, obtaining shot and keyframe information is easier than ever! Simply run the following command to get the result:
```
$ ./distribute/bin/hecate -i examples/video.mp4 --print_shot_info  --print_keyfrm_info
```

Below is the results we obtained on our dev machine (OS X 10.10 with OpenCV v3.1):
```
shots: [0:81],[84:93],[96:102],[108:270],[272:418],...,[9966:10131],[10135:10164]
keyframes: [52,85,98,128,165,208,242,259,265,273,...,10127,10141]
```
The units are frame indices (zero-based). You will notice that shot ranges are non-continuous; there are "gaps" between shots, e.g., two frames are missing between the first two shots [0:81] and [84:93]. This is normal and intentional: Hecate discards low-quality frames that aren't ideal in producing nicely looking thumbnails, animated GIFs, and video summaries. We refer to our CIKM 2016 paper for the rational behind our reason to invalidate low-quality frames.

### Thumbnail generation
Hecate uses computer vision to determine frames that are most "suitable" as video thumbnails. By suitable, we mean a frame that is the most relevant to the video content and that is the most beautiful in terms of computational aesthetics; technical details are explained in our CIKM 2016 paper.

You can generate thumbnail images using Hecate. Run the following command to generate one thumbnail image from the video.
```
$ ./distribute/bin/hecate -i examples/video.mp4 --generate_jpg --njpg 1
```
You will see the generated thumbnail image under the output directory (set as `output` by default; you can change this using the option `--out_dir YOUR_DIRECTORY`). On our dev machine we get this thumbnail image:

![alt text](https://github.com/yahoo/hecate/blob/master/examples/video_00.jpg "Hecate Thumbnail Image")

In the above example, we generated only one thumbnail image. Are you not satisfied with the thumbnail image? Hecate can generate any number of thunbmail images! Let's generate five thumbnail images.
```
$ ./distribute/bin/hecate -i examples/video.mp4 --generate_jpg --njpg 3
```

The output files are named `<video_filename>_<rank>.jpg`. The files are ranked by their quality (rank 0 means it's the best one).

### Animated GIF generation
Do you want to create animated GIFs from a video without the hassle of using manual tools? Hecate can automatically create them for you! Run the following command to create one animated GIF from the video.
```
$ ./distribute/bin/hecate -i examples/video.mp4 --generate_gif --ngif 1
```
On our dev machine, we get this animated GIF:

![alt text](https://github.com/yahoo/hecate/blob/master/examples/video_00.gif "Hecate Animated GIF")

You can, of course, create more than just one GIF by setting the paramter `--ngif N` with an appropriate number N. When there are multiple GIFs, you can also generate a "summary GIF" by concatenating them, using this command:
```
$ ./distribute/bin/hecate -i examples/video.mp4 --generate_gif --ngif 3 --generate_gifsum
```

If you'd rather want to obtain all available GIFs from the video, use the following command: 
```
$ ./distribute/bin/hecate -i examples/video.mp4 --generate_gifall
```

### Video summary generation
Last but not least, Hecate can summarize a video! Run the following command to create a video summary of length 15 seconds.
```
$ ./distribute/bin/hecate -i examples/video.mp4 --generate_mov --lmov 15
```
We included the video summary generated on our dev machine here: 
[https://github.com/yahoo/hecate/blob/master/examples/video_sum.mp4](https://github.com/yahoo/hecate/blob/master/examples/video_sum.mp4)


## Developer

Yale Song: [github](https://github.com/yalesong), [website](http://people.csail.mit.edu/yalesong)
