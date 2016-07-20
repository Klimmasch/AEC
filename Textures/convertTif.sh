#! /bin/bash
# Decompress BMP images into new folder

mkdir converted
echo -n "Converting images..."
for IMAGE in $(ls *.TIF)
do
    convert $IMAGE ./converted/${IMAGE%.*}.bmp
done
echo "done."
