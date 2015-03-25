#!/bin/bash
MAX=300
echo "clear" >> $1
echo "reset" >> $1
echo "set terminal gif animate delay 10" >> $1
#echo "set palette defined ( -1 \"black\", 0 \"blue\", 1 \"green\")" >> $1
echo "set palette rgb 23,28,3" >> $1
#echo "set pal gray" >> $1
#echo "set pm3d map" >> $1
echo "set output \""$2"\"" >> $1
echo "set isosample 40,40" >> $1
echo "set hidden3d" >> $1
echo "set xrange [0:1]" >> $1
echo "set xlabel \"X\"" >> $1
echo "set yrange [0:1]" >> $1
echo "set xlabel \"Y\"" >> $1
echo "set zrange [-0.003:0.007]" >> $1
echo "set zlabel \"Z\"" >> $1
echo "set cbrange[-0.003:0.007]" >> $1
for i in `seq 0 ${MAX}`
do
echo "splot \""$3"\" index $i using 1:2:3 with pm3d t \"waves\"" >> $1
done
