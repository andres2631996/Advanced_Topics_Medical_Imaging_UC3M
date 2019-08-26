% Example code for registration
% All the code should be executed from a folder containing both the CODE folder
% and the IMAGES folder while this code should be all in the CODE folder
close all
clear all
clc
addpath('.\CODE');
addpath('.\CODE\NIFTI_reader');
addpath('.\IMAGES');
%% Load images
file1='.\Images\Study003_MR_T1.nii';
im1 = load_nii(file1);
i1_256 = uint8(255*mat2gray(im1.img));
size_i1 = size(i1_256);
im_center1 = (size_i1+1)./2;

file2='.\Images\Study003_MR_T2.nii';
im2 = load_nii(file2);
i2_256=uint8(255*mat2gray(im2.img));
size_i2=size(im2.img);
im_center2 = (size_i2+1)./2;

%% DEFINE TRANSFORMATION PARAMETERS
% The following line specifies the transformation that will be applied to
% image2. The registration algorithm will try to find the parameters that
% we have used to transform the algorithm.
% The first 3 parameters are rotations in degrees, the next 3 are
% translations and then three scalings (all 1.0)
Tr_Mat = makeTransf_3D_center(20,40,10,-10,-10,-10,1,1,1,im_center2); % MODIFY THIS LINE TO TEST DIFFERENT TRANSFORMATIONS
i2_Transf = transform_image_3D(Tr_Mat,single(i2_256), 'linear');
i2_Transf_256 = uint8(255*mat2gray(i2_Transf));

%% OPTIMIZATION
% The following lines will perform the registration step. The following 
% parameteres will initialize the optimization. We don't have an initial 
% guess of the transformation, so they are zero (no rotations and no translations)
params_ini = [0 0 0 0 0 0]; % initial parameters
samp = [2 2 2];

fi=fopen('AlgorithmValues.txt','w+');
diary('AlgorithmValues.txt');
diary on
tic
[p_final,f_value] = spm_powell(params_ini,eye(6),0.001,'cost_f',i1_256,i2_Transf_256, samp); % powell optimizer
toc



%% DISPLAY REGISTRATION RESULT
disp(' ')
disp(['Transformation (FINAL): ', num2str(p_final')])
disp(['Minus H (FINAL): ', num2str(f_value)])
diary off
%% APPLY REGISTRATION RESULT
% The following lines will apply the transformation found in the registration step
Tr_Mat2 = makeTransf_3D_center(p_final(1),p_final(2),p_final(3),p_final(4),p_final(5),p_final(6),1,1,1,im_center2);
i2_Transf_TransfBack = transform_image_3D(inv(Tr_Mat2),single(i2_Transf_256), 'linear');
i2_Transf_TransfBack_256 = uint8(255*mat2gray(i2_Transf_TransfBack));

% Display
figure
subplot(2,2,1), imshow(squeeze(i1_256(:,:,floor(im_center1(3)))),[])
subplot(2,2,1), title('Image 1 (T1)')
subplot(2,2,2), imshow(squeeze(i2_Transf_256(:,:,floor(im_center2(3)))),[])
subplot(2,2,2), title('Image 2 (T2): Not registered')
subplot(2,2,3:4), imshow(squeeze(i2_Transf_TransfBack_256(:,:,floor(im_center2(3)))),[])
subplot(2,2,3:4), title('Image 2 registered')

% Saving all the displayed parameters on screen in a .txt file and making a
% plot with them
rx=zeros;
ry=zeros;
rz=zeros;
tx=zeros;
ty=zeros;
tz=zeros;
cost_f=zeros;
cont_param=0;
cont_f=0;
while ~feof(fi)
    line=fgetl(fi);
    if isempty(strfind(line,'Transformation:'))==0 || isempty(strfind(line,'Transformation (FINAL):'))==0 || isempty(strfind(line,'Minus H:'))==0 || isempty(strfind(line,'Minus H (FINAL):'))==0
        if strcmp(line(1:8),'Minus H:')==1
            cont_f=cont_f+1;
            cost_f(cont_f)=str2num(line(10:end));
        elseif strcmp(line(1:15),'Transformation:')==1
            cont_param=cont_param+1;
            x=str2num(line(17:end));
            rx(cont_param)=x(1);
            ry(cont_param)=x(2);
            rz(cont_param)=x(3);
            tx(cont_param)=x(4);
            ty(cont_param)=x(5);
            tz(cont_param)=x(6);
        elseif strcmp(line(1:16),'Minus H (FINAL):')==1
            cost_f(end)=str2num(line(18:end));
        elseif strcmp(line(1:23),'Transformation (FINAL):')==1
            x=str2num(line(25:end));
            rx(end)=x(1);
            ry(end)=x(2);
            rz(end)=x(3);
            tx(end)=x(4);
            ty(end)=x(5);
            tz(end)=x(6);
        end
    end
end
figure
subplot(2,3,1),plot(1:cont_param,rx),title('Rotation in X'), xlabel('Iteration number'), ylabel('Degrees')
subplot(2,3,2),plot(1:cont_param,ry),title('Rotation in Y'), xlabel('Iteration number'), ylabel('Degrees')
subplot(2,3,3),plot(1:cont_param,rz),title('Rotation in Z'), xlabel('Iteration number'), ylabel('Degrees')
subplot(2,3,4),plot(1:cont_param,tx),title('Translation in X'), xlabel('Iteration number'), ylabel('Voxels')
subplot(2,3,5),plot(1:cont_param,ty),title('Translation in Y'), xlabel('Iteration number'), ylabel('Voxels')
subplot(2,3,6),plot(1:cont_param,tz),title('Translation in Z'), xlabel('Iteration number'), ylabel('Voxels')
figure, plot(1:cont_f,-cost_f),title('Variation of Mutual Information in the registration  algorithm'),xlabel('Iteration number'),ylabel('Mutual Information (MI)')
fclose(fi);