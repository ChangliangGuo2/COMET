%--------------------------------------------------------------------------
clear all;
close all;
clc;
%--------------------------------------------------------------------------
%%Load C traces values obtained from CaImAn pipeline
%%C traces refer to the cellnumber X framenumber; for example if you get 1000 cells from CaImAn
%%and have 10000 frames the C traces matrix should be 1000X10000
Folder     = 'D:\Cscope Data\PP367\Session1_12_51_06\Processed\';
load([Folder 'C_Raw.mat']);
C          = single(C_Raw);
%--------------------------------------------------------------------------
%%Load Registed cells from CaImAn, if you get 1000 cells and the frame size is 600X600, then 
%%A neuron is 600X600X1000
load([Folder 'A_T.mat']);
A_Stack = single(A);
%--------------------------------------------------------------------------
%%Load the index of the good cells reconstructed
load([Folder 'A_neuron_good_idx.mat']);
A_neuron_good = foo;
A_neuron_good = single(A_neuron_good)+1;
%--------------------------------------------------------------------------
%%Only choose the good cells in A and C traces in 1000 cells; for example, 
%%there are 800 good cells in 1000 cells
A_Stack_good = A_Stack(:,:,A_neuron_good);
C_good       = C;
%--------------------------------------------------------------------------
%%Reshape the 3D matrix A neurons from 3D into 2D; for example, if good A neuron
%%is 600X600X800 then the 2D is 360000X800
A_Stack_2D = reshape( A_Stack_good, size(A_Stack_good,1)*size(A_Stack_good,2),size(A_Stack_good,3) );
%--------------------------------------------------------------------------
%% Show colored cells
K    = size(C_good,1);
temp = prism; % prism is the a 256X3 color map
% temp = bsxfun(@times, temp, 1./sum(temp,2));
col  = temp(randi(64, K,1), :); % randi(64,K,1) generats KX1 matrix pseudorandom integers between 1 and 64.
size(diag(col(:,1))*C_good)
for m=1:3
    Y_mixed(:, :, m) = A_Stack_2D* (diag(col(:,m))*C_good);
    m
end
figure(100)
imshow(uint8(C_good ))
%--------------------------------------------------------------------------
%%Save colored cells with frames into tif stack file
Y_mixed_4D = reshape( Y_mixed, size(A_Stack_good,1), size(A_Stack_good,2),size(C_good,2) ,3 );
for Fr_index = 1:size(C_good,2)
    Y_mixed_Stack(:,:,:) = reshape( Y_mixed_4D(:,:,Fr_index,:) , size(A_Stack_good,1), size(A_Stack_good,2), 3);
    Y_mixed_Stack = Y_mixed_Stack*(255/max(Y_mixed_Stack(:)));
    imwrite(uint8(Y_mixed_Stack),[Folder 'Y_mixed_Stack_test_nomask.tif'],'WriteMode','append')
    disp( [num2str(Fr_index) '/' num2str(size(C,2))]);
end
%--------------------------------------------------------------------------