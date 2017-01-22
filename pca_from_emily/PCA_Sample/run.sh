#! /bin/sh

export PATH=$PATH:/home/naoki/Desktop/pcasift-0.91p/
find . -name "*.pgm" -exec /home/naoki/Desktop/pcasift-0.91p/pcasift /home/naoki/Desktop/pcasift-0.91p/gpcavects.txt {} {}.txt \;
