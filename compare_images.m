function [diff,H,i2_256]=compare_images(i1,i2,Tr_Mat,costF)

% This function takes two images (i1 and i2) and compares i1 to a transformed version
% of i2 applying Tr_M transformation Matrix

i2_Tr = transform_image_3D(Tr_Mat,single(i2.img), 'linear');

i1_256=uint8(255*mat2gray(i1.img));
i2_256=uint8(255*mat2gray(i2_Tr));


switch costF
    case 'Corr'
        % Compute correlation
        diff =corr2(i1_256(:),i2_256(:));
    case 'Abs_Diff'
        % Compute absolute differences
        abs_diff=abs(i1_256-i2_256);
        diff =sum(abs_diff(:))/(size(abs_diff,1)*size(abs_diff,2)*size(abs_diff,3));
    case 'MI'
        % Compute joint histogram
        H = spm_hist2(i1_256,i2_256, eye(4) ,[1 1 1]);
        
        % Compute cost function from histogram
        H  = H+eps;
        H  = H/sum(H(:)); % normalization
        s1 = sum(H,1); % marginal probability image 1
        s2 = sum(H,2); % marginal probability image 2
        
        %  Compute mutual Information:
        entropy1=-sum(s1.*log2(s1));
        entropy2=-sum(s2.*log2(s2));
        joint_entropy=-sum(sum(H.*log2(H)));
		diff  =entropy1+entropy2-joint_entropy;
        
        
    otherwise
        diff=0;
end

end
