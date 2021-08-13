% Little script that enables zooming in and out of the images to
% simulate aniseikonia without changing the size of the images.
function img_out = zoomImgKeepSize(img_in, factor)

      if factor == 0
          img_out = [];
          return;
      end
      if factor < 0
          img_in = rot90(img_in, 2);
          factor = 0 - factor;
      end
      if factor == 1
          img_out = img_in;
          return;
      end
      xSize = size(img_in, 1);
      ySize = size(img_in, 2);
      xCrop = ceil(xSize / 2 * abs(factor - 1));
      yCrop = ceil(ySize / 2 * abs(factor - 1));
      zoomPicture = imresize(img_in, factor);
      if factor > 1
          img_out = zoomPicture( 1+xCrop:end-xCrop, 1+yCrop:end-yCrop );
      else
          img_out = zeros(xSize, ySize);
          img_out( 1+xCrop:end-xCrop, 1+yCrop:end-yCrop ) = zoomPicture;
      end

end
