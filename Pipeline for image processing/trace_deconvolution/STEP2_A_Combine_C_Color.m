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
load([Folder 'A_T_clean.mat']);
load([Folder 'C_Raw_CaAct.mat']);
C_good_list = ~ismember(1:size(C_Raw,1),A_bad_list);
C_Raw_clean_CaAct = C_Raw( C_good_list, : );
%--------------------------------------------------------------------------
%%Load the index of the good cells reconstructed
%--------------------------------------------------------------------------
%%Only choose the good cells in A and C traces in 3675 cells; for example, 
%%there are 1772 good cells in 3675 cells
A_good = A_good_3D;
A_x = size(A_good,1);
A_y = size(A_good,2);
A_z = size(A_good,3);

A_good_2D = reshape( A_good , A_x*A_y , A_z );
A_good_2D = single(A_good_2D);
C = single(C_Raw_clean_CaAct);
K = size(C,1);

temp = prism;
color = temp(randi(256,K,1),:);

FrameNo = 1000;
Framleft = 1000;%size(C,2);

for m=1:3
    C_colored(:, :, m) = (diag(color(:,m))*C);
end

save([Folder 'C-raw_C-color.mat'],'C_Raw_clean', 'C_colored','-v7.3');

FramStart = 1;
while (Framleft>0)

    C_sub = C(:,FramStart:FramStart-1 + min(FrameNo, Framleft));
    Y_mixed = zeros( size(A_good_2D,1) , min(FrameNo, Framleft) ,3 );
    for m=1:3
        Y_mixed(:,:,m)=A_good_2D*(diag(color(:,m))*C_sub);       
    end
    figure(100);
    imshow(uint8(C_sub));
    Y_mixed_4D = reshape( Y_mixed , A_x,A_y,size(C_sub,2),3 );

    for Fram_idx = 1:size(C_sub,2)
        Y_mixed_temp = reshape( Y_mixed_4D(:,:,Fram_idx,:),A_x,A_y,3 );
        Y_mixed_temp = Y_mixed_temp*(255/max(Y_mixed_temp(:)));
        imwrite(uint8(Y_mixed_temp),[Folder 'ColorCells_3D_' num2str(1) '-' num2str(13085) '.tif'],'WriteMode','append')
        
        disp([num2str(Fram_idx) '/' num2str(size(C_sub,2))]);
    end
    
    FramStart = FramStart + FrameNo;
    Framleft = Framleft-FrameNo;
    
    disp([num2str(Framleft) '/' num2str(size(C,2))]);

end

C_Raw_Norm = normalize(C_Raw_clean,2,'range');
    for grp=1:17
         figure(200+grp);
    for idx=1:50
       
        plot(C_Raw_Norm(idx+(grp-1)*50,:)+idx);
        hold on
    end
    end

R_v=[]
P_v=[]
for idx=1:size(C_Raw_clean,1)

    [R,P] = corrcoef(C_Raw_clean(1,:),C_Raw_clean(idx,:));
    R_v(idx)=R(1,2);
    P_v(idx)=P(1,2);
    idx

end
figure(1000)
plot(R_v)
hold on;
plot(P_v)
C_sim_list = find(R_v>=0.45)
figure(1001)
    for idx=1:length(C_sim_list)
       
        plot(C_Raw_Norm(idx,:)+idx);
        hold on
    end


A_similar_C = normalize(A_good_2D(:,C_sim_list),1,'range');
  A_color_similar_C = (A_similar_C.*(A_similar_C>0.4))*color(C_sim_list,:);
A_color_2D_similar_C = reshape( A_color_similar_C , A_x,A_y,3 )*1;




A_good_2D_norm = normalize(A_good_2D,1,'range');
  A_color = (A_good_2D_norm.*(A_good_2D_norm>0.4))*color;
A_color_2D = reshape( A_color , A_x,A_y,3 )*1;