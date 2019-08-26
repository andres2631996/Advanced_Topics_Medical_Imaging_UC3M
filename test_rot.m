
function [diff_array1,diff_array2,diff_array3]=test_rot(i1,i2)

size_i1=size(i1.img);
im_center1 = (size_i1+1)./2;

% Definition
angle=-35:5:35;
diff_array1=zeros(1,length(angle));
diff_array2=zeros(1,length(angle));
diff_array3=zeros(1,length(angle));
% Rotation + Cost function value computation
for i = 1:length(angle)
    Tr = makeTransf_3D_center(0,0,angle(i),0,0,0,1,1,1,im_center1);
    diff_array1(i)=compare_images(i1,i2,Tr,'Abs_Diff');
    [diff_array2(i),H,i2_256]=compare_images(i1,i2,Tr,'MI');
    diff_array3(i)=compare_images(i1,i2,Tr,'Corr');
    h=figure;
    title('Rotations');
    subplot(3,3,1),imshow(squeeze(i1.img(:,:,round(size(i1.img,3)/2))),[]);
    subplot(3,3,2),imshow(squeeze(i1.img(:,round(size(i1.img,2)/2),:)),[]);
    subplot(3,3,3),imshow(squeeze(i1.img(round(size(i1.img,1)/2),:,:)),[]);
    subplot(3,3,4),imshow(squeeze(i2_256(:,:,round(size(i2_256,3)/2))),[]);
    subplot(3,3,5),imshow(squeeze(i2_256(:,round(size(i2_256,2)/2),:)),[]);
    subplot(3,3,6),imshow(squeeze(i2_256(round(size(i2_256,1)/2),:,:)),[]);
    subplot(3,3,8),imshow(log(H),[]),title('Logarithmic joint histogram');
    pause(2);
    close(h);
end

% Generate plots
figure
plot(angle, diff_array1, 'b');
xlabel('Angle (�)')
ylabel('Sum of Absolute Differences')

figure
plot(angle, diff_array2, 'g');
xlabel('Angle (�)')
ylabel('Mutual Information')

figure
plot(angle, diff_array3, 'r');
xlabel('Angle (�)')
ylabel('Correlation')

end
