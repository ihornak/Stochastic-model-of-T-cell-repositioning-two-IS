#!/bin/bash

rm povrayFiles/*.pov
rm pictures/*.png

echo "Treatment of text files and preparation of povray files"
python Make_Povray_files.py


echo "Making of pictures"
for i in  povrayFiles/*.pov
do
    if test -f "$i" 
    then
       povray -D Height=562 Width=750 "$i"
    fi
done



echo "Moving pictures"
mv povrayFiles/*.png  pictures/


ffmpeg -framerate 2/1 -start_number 0 -i pictures/microtubule_%1d.png -c:v libx264 -r 30 -pix_fmt yuv420p two_IS.mp4

