#!/bin/bash

if [ $# -lt 1]; then
	echo usage $0 file_base
	exit 1;

fi

ffmpeg -framerate 30 -i $1_%d.png -c:v libx264 -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" -profile:v high -crf 20 -pix_fmt yuv420p $1.mp4
