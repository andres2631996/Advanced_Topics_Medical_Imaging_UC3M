% Example code
% All the code should be executed from a folder containing both the CODE folder
% and the IMAGES folder while this code should be all in the CODE folder
close all
clear all
clc
addpath('.\CODE');
addpath('.\CODE\NIFTI_reader');
addpath('.\IMAGES');
%% Define image names
file1='.\IMAGES\MR.nii';
file2='.\IMAGES\PET.nii';

%% Load images using Nifty reader
im1=load_nii(file1);
im2=load_nii(file2);

%% Test rotation/translation
[d1rot,d2rot,d3rot]=test_rot(im1,im2);
[d1trans,d2trans,d3trans]=test_trans(im1,im2);
