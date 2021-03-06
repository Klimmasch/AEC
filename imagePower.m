%%% Small script for analyzing image properties

close all;

figure;
%% good ones:
img = imread('Textures/vanHateren/17.bmp');
% img = imread('Textures/vanHateren/40.bmp');
% img = imread('Textures/vanHateren/85.bmp');
% img = imread('Textures/vanHateren/cutOffScale_corrected_2219.bmp');
% img = imread('Textures/vanHateren/cutOffScale_corrected_2473.bmp');

%% bad ones:
% img = imread('Textures/vanHateren/suboptimal/6.bmp');
% img = imread('Textures/vanHateren/suboptimal/cutOffScale_corrected_2894.bmp');
% img = imread('Textures/vanHateren/suboptimal/47.bmp');
% img = imread('/home/lelais/Documents/MATLAB/aec/Textures/mcgillManmade/pippin_win07_017.bmp');
% imagesc(img);

% shift 0-frequency coefficient from left upper corner into center
imgs = fftshift(img(:, :, 2));

% fourier spectrum
im_fft = fft2(imgs);

% power
im_Pfft = abs(fftshift(im_fft)) .^ 2;

% to brighten the display use log of power, avoid log of zeros
im_logPfft = log(im_Pfft + eps);

% scaling for better visualization
im_sclogPfft = im_logPfft * 100;

figure;
imagesc(im_logPfft);
% colormap(gray);
colorbar;
title('magnitude spectrum');

figure;
imagesc(angle(im_fft));
colormap(gray);
colorbar;
title('phase spectrum');

figure;
histogram(im_logPfft, 100, 'Normalization', 'pdf');
xlabel('Power');
ylabel('Probability');

% figure;
% histogram(angle(im_fft), 100, 'Normalization', 'pdf');
% xlabel('Phase');
% ylabel('Probability');

% Frequency Domain Filtering of image
hz = fspecial('sobel');
PQ = paddedsize(size(img));
HZ = fft2(double(hz), PQ(1), PQ(2));
F = fft2(double(img), PQ(1), PQ(2));
FDF = HZ .* F(: , :, 1);
fdf = ifft2(FDF);
fdf = fdf(1 : size(img, 1), 1 : size(img, 2));

% use absolute function to get rid of negative frequencies
figure;
imshow(abs(fdf), []);
colorbar;


% threshold into a binary image
figure;
imshow(abs(fdf) > 0.1 * abs(max(fdf(:))));
colorbar;
