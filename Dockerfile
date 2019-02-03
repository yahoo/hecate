#FROM valian/docker-python-opencv-ffmpeg
FROM m03geek/ffmpeg-opencv:stretch

RUN apt-get update && apt-get -yqq install wget unzip make g++
RUN wget -q https://github.com/rllin-fathom/hecate/archive/master.zip -O hecate.zip\
    && unzip -qq hecate.zip \
    && mv hecate-master hecate \
    && cd hecate \
    && make all \
    && make distribute
