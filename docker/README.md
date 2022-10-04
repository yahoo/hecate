# HECATE - DOCKER INSTALLATION


## Installation

First you need to install docker on your local computer, see following [tutorial](https://docs.docker.com/install/linux/docker-ce/ubuntu/#set-up-the-repository). Note, for running the docker properly you have be logged as superuser otherwise you will face many partial issues which sometimes does not make much sense.

Check your installation using
```
docker -v
```

Once you install the docker correctly, follow the instruction below:
```
$ git clone https://github.com/yahoo/hecate.git
$ cd hecate/docker
$ ./build.sh
```

Now the hecate CLI should open up for you or else you may need to:
```
docker exec -it hecate bash
```

Once you've successfully compiled hecate, it will generate a binary executable under `distribute/bin/`. Run the following command to check if everything works properly:
```
$ cd hecate
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


## Docker Developer

Sayantan Mandal: [github](https://github.com/smandal047)

## 