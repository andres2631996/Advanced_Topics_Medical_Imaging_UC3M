function[dice]=dice_coefficient(I_seg,I_seg_diff)
% Computes the Dice coefficient between two segmented images, looking for
% non-zero values in both segmented images, assuming that zero values make
% part of the background
ind_1=find(I_seg);
ind_2=find(I_seg_diff);
intersection=intersect(ind_1,ind_2);
dice=2*length(intersection)/(length(ind_1)+length(ind_2));
