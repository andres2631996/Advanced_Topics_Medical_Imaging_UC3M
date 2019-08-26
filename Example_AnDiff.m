%% Load image
I_orig  = imread('MR_axial_noise.tif');

%% Define parameters
im=I_orig;
num_iter = 50;
delta_t = 1/5;
kappa = 15;
option = 1;

%% Apply anisotropic diffusion function
tic
diff_im = anisodiff2D(im, num_iter, delta_t, kappa, option);
toc
%% Display images
figure, imshow(I_orig,[]);
figure, imshow(diff_im,[]);
