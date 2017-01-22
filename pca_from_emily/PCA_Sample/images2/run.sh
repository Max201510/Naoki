#! /bin/sh

export PATH=$PATH:/home/naoki/Desktop/pcasift-0.91p/images2
find . -name "*.pgm" -exec ../pcasift ../gpcavects.txt {} {}.pkeys \;


