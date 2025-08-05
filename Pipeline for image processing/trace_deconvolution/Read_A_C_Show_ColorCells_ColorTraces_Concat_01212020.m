close all;
clear all;
clc;

% Folder = '/home/chang/caiman_data/example_movies/Caiman_1/';
Folder = 'M:\2020_12_22 @ wirelessbehav\16_35_34\experiment\12-22-2020\SD2\Caiman_1\';

load([Folder 'A_neuron_1_3d_arry.mat']);
A = foo;


load([Folder 'Video1-14_8bit_els__d1_608_d2_608_d3_1_order_F_frames_13085_C.mat']);
C = foo;

clear foo;

A_x = size(A,1);
A_y = size(A,2);
A_z = size(A,3);

A_2D = reshape( A , A_x*A_y , A_z );
A_2D = single(A_2D);
clear A;
C = single(C);
K = size(C,1);

temp = prism;
color = temp(randi(256,K,1),:);

FrameNo = 2000;
Framleft = size(C,2);

    for m=1:3
        C_colored(:, :, m) = (diag(color(:,m))*C);
    end

save([Folder 'C_colored.mat'],'C_colored','-v7.3');

FramStart = 1;
while (Framleft>0)

    C_sub = C(:,FramStart:FramStart-1 + min(FrameNo, Framleft));
    Y_mixed = zeros( size(A_2D,1) , min(FrameNo, Framleft) ,3 );
%     Y_mixed = 0;
    for m=1:3
        Y_mixed(:,:,m)=A_2D*(diag(color(:,m))*C_sub);
       
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