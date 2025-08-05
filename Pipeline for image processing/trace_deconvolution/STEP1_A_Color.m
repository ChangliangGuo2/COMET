%--------------------------------------------------------------------------
%-----------------This function is used to concatenate the A matrix 
clear all;
close all;
clc;
%-----------------By Changliang Guo at 10032022
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%----
Folder       = 'D:\Cscope Data\Mice_New_LinearTrack\PP367_Ses1_11_09_30\Processed\';
A_matFiles   = dir([Folder 'A_neuron_1_3d_arry*.mat']);
A            = [];
for idx      = 1:numel(A_matFiles)
    A_names  = A_matFiles(idx).name;
    A_var    = struct2cell(load([Folder A_names]));
    A_var    = A_var{1};
    A        = cat(3,A,A_var);
    disp([num2str(idx) '/' num2str(numel(A_matFiles))])
end
%--------------------------------------------------------------------------
%%Load the index of the good cells reconstructed
load([Folder 'A_neuron_good_idx.mat']);
A_listgood = foo;
A_listgood = single(A_listgood)+1;
%--------------------------------------------------------------------------
A_good = A(:,:,A_listgood);
A_x = size(A_good,1);
A_y = size(A_good,2);
A_z = size(A_good,3);

A_good_2D = reshape( A_good , A_x*A_y , A_z );
A_good_2D = single(A_good_2D);
K = size(A_good,3);

temp = prism;
color = temp(randi(256,K,1),:);



A_good_2D_norm = normalize(A_good_2D,1,'range');
A_color = (A_good_2D_norm.*(A_good_2D_norm>0.4))*color;
A_color_2D = reshape( A_color , A_x,A_y,3 )*1;

figure(100)
imshow(A_color_2D)                                                                                                                                                                                                                                                                                                                                                                                                                                      




for idx=1:size(A_good_2D,2)
    figure(200)
    imshow(A_good(:,:,idx)*4)
    idx
pause()
end

A_bad_list = [1 8 12 16 22 24 35 52 103 106 118 122 124 129 131 191 199 200 201 206 209 213 214 216 276 282 310 325 441 449 461 476 479 487 492 493 529 567 568 571 575 585 587 590 591 601 607 618 633 652 653 694 732 735 736 739 740 742 743 757 780 814 827 831 832 837 846 860 861 864 873 903 904 910 927 929 930 943 ]
A_good_list = ~ismember(1:size(A_good,3),A_bad_list);
A_good_2D_norm = normalize(A_good_2D(:,A_good_list),1,'range');
K = size(A_good_2D_norm,2);
temp = prism;
color = temp(randi(256,K,1),:);
Threshold=0.4;
A_color = (A_good_2D_norm.*(A_good_2D_norm>Threshold))*color;
A_color_2D = reshape( A_color , A_x,A_y,3 )*1;

figure(300)
imshow(uint8(A_color_2D*255))

imwrite(uint8(A_color_2D*255),[Folder 'A_Color_Total_Rand_Trd' num2str(Threshold) '.tif']);

A_good_3D = reshape( A_good_2D_norm , A_x,A_y,K )*1;

% for idx=1:size(A_good_3D,3)
%     figure(200)
%     imshow(A_good_3D(:,:,idx)*4)
%     idx
% pause()
% end

save([Folder 'A_bad_list_clean.mat'],'A_bad_list','-v7.3');

%------------Load C traces for showing the dynamics of the activity
load([Folder 'C_Raw.mat']);
%-----------------------------
C_Raw_clean = C_Raw(A_good_list,:);
save([Folder 'C_Raw_clean.mat'],'C_Raw_clean','-v7.3');

save([Folder 'A_T_clean.mat'],'A_good_3D','A_bad_list','C_Raw_clean','-v7.3');