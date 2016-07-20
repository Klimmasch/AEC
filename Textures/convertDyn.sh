#! /bin/bash
# Converts images of various formats into BMP images
# into a newly created folder and optionally cuts them
# into a qudratic resultion and handling offsets dynamically

# test number of passed command line arguments
if [[ $# -lt 2 ]]
then
    echo "Error: Not enough imput arguments provided!"
    echo "Usage: $(basename $0) RAW_IMAGE_DIR CONVERTED_IMAGE_DIR [CUT]"
    echo "CUT=[1]: whether to cut the converted image into a quadratic resolution"
    echo ""
    exit 1
fi

CUR_DIR="$(echo $(pwd))"
RAW_IMAGE_DIR="$1"
CONVERTED_IMAGE_DIR="$2"

mkdir $CONVERTED_IMAGE_DIR > /dev/null 2>&1

if [[ $# -eq 2 ]]
then
    # just convert
    echo -n "Converting images..."
    for IMAGE in $(ls $RAW_IMAGE_DIR)
    do
        convert "$RAW_IMAGE_DIR/$IMAGE" "$CONVERTED_IMAGE_DIR/${IMAGE%.*}.bmp"
    done
else
    # convert & cut
    echo -n "Converting and cutting images..."
    for IMAGE in $(ls $RAW_IMAGE_DIR)
    do
        NEW_IMG="${IMAGE%.*}.bmp"
        convert "$RAW_IMAGE_DIR/$IMAGE" "$CONVERTED_IMAGE_DIR/$NEW_IMG"

        WIDTH=$(convert $CONVERTED_IMAGE_DIR/$NEW_IMG -print "%w" /dev/null)
        HEIGHT=$(convert $CONVERTED_IMAGE_DIR/$NEW_IMG -print "%h" /dev/null)
        if [[ $HEIGHT -gt $WIDTH ]]
        then
            convert "$CONVERTED_IMAGE_DIR/$NEW_IMG" -crop "$WIDTH""x""$WIDTH+0+0" "$CONVERTED_IMAGE_DIR/$NEW_IMG"
        else
            OFFSET=$(($(($WIDTH - $HEIGHT)) / 2))
            convert "$CONVERTED_IMAGE_DIR/$NEW_IMG" -crop "$HEIGHT""x""$HEIGHT+$OFFSET+0" "$CONVERTED_IMAGE_DIR/$NEW_IMG"
        fi
    done
fi
echo "done."

exit 0
