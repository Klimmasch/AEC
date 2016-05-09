#! /bin/bash
# Decompress BMP images into new folder

mkdir decompressed
echo -n "Decompressing images..."
for IMAGE in $(ls *.bmp)
do
    convert $IMAGE +matte +compress ./decompressed/$IMAGE
done
echo "done."
