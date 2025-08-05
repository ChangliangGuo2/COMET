%--------------------------------------------------------------------------
%-----------------This function is used to concatenate the A matrix 
clear all;
close all;
clc;
%-----------------By Changliang Guo at 10032022
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%----
Folder       = 'D:\Cscope Data\PP367\Session1_12_51_06\Processed\';
load([Folder 'A_T_clean.mat']);

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
C = single(C_Raw_clean);
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

C_Raw_Norm = normalize(C_Raw,2,'range');
    figure(200);
    for idx=1:50
        plot(C_Raw_Norm(idx,:)+idx);
        hold on
    end

A_good_2D_norm = normalize(A_good_2D,1,'range');
A_color = (A_good_2D_norm.*(A_good_2D_norm>0.4))*color;
A_color_2D = reshape( A_color , A_x,A_y,3 )*1;