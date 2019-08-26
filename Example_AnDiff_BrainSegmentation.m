clc, close all, clear all
%% Load image
I = imread('Brain.tif');

%% Add noise
I_noisy = imnoise(I,'gaussian',0,0.015);

%% Define parameters
im=I_noisy;
num_iter = 100;
delta_t = 1/5;
kappa = 15;
option = 1;

%% Apply anisotropic diffusion function
diff_im = anisodiff2D(im, num_iter, delta_t, kappa, option);

%% Display images
figure(1)
subplot(1,3,1)
imshow(I,[0,255]);
title('Original')
subplot(1,3,2)
imshow(I_noisy, [0,255]);
title('Noisy')
subplot(1,3,3)
imshow(diff_im,[0,255]);
title('Diffused')

%% Segmentation
[M,N]=size(I);
I = double(I);
I_seg = zeros(M,N);
pos = find(I>=200);
I_seg(pos) = 1;

I_seg_diff = zeros(M,N);
pos = find(diff_im>=200);
I_seg_diff(pos) = 1;

figure(2)
subplot(1,2,1)
imshow(I_seg,[])
title('Original')
subplot(1,2,2)
imshow(I_seg_diff,[])
title('Diffused')

dice_coefficient(I_seg,I_seg_diff)
 
