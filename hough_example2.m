close('all')


I  = imread('test_line2_small.tif');
% This time we do not rotate the image
rotI = rgb2gray(I(:,:,1:3));

fig1 = imshow(rotI);

%Find the edges in the image using the edge function.
BW = edge(rotI,'canny');
figure, imshow(BW);

%Compute the Hough transform of the image using the hough function or your own function HoughTransform_basic.
[H1,theta1,rho1] = hough(BW,'RhoResolution',5.0,'ThetaResolution',5.0);

%Display the transform using the imshow function.
figure, imshow(imadjust(mat2gray(H1)),[],'XData',theta1,'YData',rho1,...
        'InitialMagnification','fit');
xlabel('\theta (degrees)'), ylabel('\rho');
title('Hough Transform Matlab');
axis on, axis normal, hold on;
colormap(hot)

%Compute the Hough transform of the image using your %own function HoughTransform_basic.
[rho2,theta2,H2] = HoughTransform_fast(BW,5.0,5.0);
%Display the transform using the imshow function.
figure, imshow(imadjust(mat2gray(H2)),[],'XData',theta2,'YData',rho2,...
        'InitialMagnification','fit');
xlabel('\theta (degrees)'), ylabel('\rho');
title('My Hough Transform');
axis on, axis normal, hold on;
colormap(hot)