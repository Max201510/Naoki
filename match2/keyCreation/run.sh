#! /bin/sh

#export PATH=$PATH:/home/naoki/Desktop/keyCreation
export PATH=$PATH:/home/max/Naoki/match2/keyCreation
find . -name "*.pgm" -exec pcasift gpcavects.txt {} {}.pkeys \;


