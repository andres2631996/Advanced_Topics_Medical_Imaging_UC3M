%% Load image
I_orig  = imread('MR_sagittal_noise.tif');

%% Define parameters
im=I_orig;
num_iter = 50;
delta_t = 1/5;
kappa = 30;
option = 1;

%% Apply anisotropic diffusion function
diff_im = anisodiff2D(im, num_iter, delta_t, kappa, option);
diff_im_linear=anisodiff2D_linear(im, num_iter, delta_t, kappa, option);
%% Display images
figure, imshow(I_orig,[]);
figure, imshow(diff_im,[]);
figure, imshow(diff_im_linear,[]);