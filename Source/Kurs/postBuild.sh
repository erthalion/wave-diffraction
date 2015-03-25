#!/bin/bash
bash create_animation.sh $1.gpi $1.gif $1.dat
gnuplot $1'.gpi'
bash create_avi.sh $1.gif $1.avi
