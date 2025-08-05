function [FR_4D_gau, info_Shuffle,infoP] = Info_Gau_ClG_Shuffle1st_Placefield_Fig7(posx, posy, Speed_r_th_idx, Dt_spk,occThresh, bin_edges, spks, sigma, numShuffles)

binSize                           = bin_edges(2)-bin_edges(1); %% bin size is 2cm
Sigma                             = sigma;
numShuffles                       = numShuffles+1;
ns                                = size(spks,1); %% number of cells                 = numShuffles+1;
info                              = nan(numShuffles,ns);%% number of shuffles X number of cells

posx                              = posx(Speed_r_th_idx);
posy                              = posy(Speed_r_th_idx);
NoFrames                          = length(posx);%% changed to the number of the frame after speed threshold, not the total number of frames.

shift                             = floor( rand(numShuffles,1) *NoFrames);
shift(1,1)                        = 0;
indexp                            = 1:1:NoFrames;
Dt_spk                            = Dt_spk'; 
dt_m                              = median(Dt_spk); 
%--------------------------------------------------------------------------
%% Circular mask for the shape of open field after zeropadding
CirX                              = bin_edges(1):0.1:bin_edges(end);
CirY                              = CirX;
[Cols Rows]                       = meshgrid(CirX, CirY);
centerX                           = 0;
centerY                           = 0;
Radius                            = bin_edges(end);
Cols_1D                           = reshape(Cols,[],1);
Rows_1D                           = reshape(Rows,[],1);
Mask_Cir_2d                       = (Cols - centerY).^2 + (Rows - centerX).^2 <= Radius.^2;
Mask_Cir                          = (Cols_1D - centerY).^2 + (Rows_1D - centerX).^2 <= Radius.^2;
Cols_P                            = Cols(Mask_Cir);
Rows_P                            = Rows(Mask_Cir);

[~, ~, ~, xbin_cir, ybin_cir]     = histcounts2(Cols_P, Rows_P, bin_edges, bin_edges);
bin2_cir_1D                       = (xbin_cir-1)*(size(bin_edges,2)-1)+ybin_cir;
nbins_cir                         = accumarray(bin2_cir_1D,1);
nbins_cir_2d                      = zeros((length(bin_edges)-1)*(length(bin_edges)-1),1); %% Bin No. X 1
nbins_cir_2d(1:length(nbins_cir)) = nbins_cir;
nbins_cir_2d_2d                   = reshape(nbins_cir_2d,[(length(bin_edges)-1),(length(bin_edges)-1)]);
nbins_cir_2d_2d                   = double(nbins_cir_2d_2d>0);
nbins_cir_2d_2d_pad               = padarray(nbins_cir_2d_2d,[5 5],0,'both');

posx2                             = zeros(NoFrames,1);
posy2                             = zeros(NoFrames,1);
%--------------------------------------------------------------------------
%----------Occupancy time
occ                               = zeros((length(bin_edges)-1)*(length(bin_edges)-1),1); %% Bin No. X 1
[~, ~, ~, xbin2, ybin2]           = histcounts2(posx, posy, bin_edges, bin_edges);
% ---------------------------------------------------
bin2_1D                           = (xbin2-1)*(size(bin_edges,2)-1)+ybin2;% Change the 2D bin indexes to 1D indexes for quick processing
% ---------------------------------------------------
temp                              = accumarray(bin2_1D,dt_m);
occ(1:length(temp))               = temp; %% occupantion time in each bin!!!
occ(isnan(occ))                   = 0;

occ_2D                            = reshape(occ',[(length(bin_edges)-1),(length(bin_edges)-1)]);
occ_2D_pad                        = padarray(occ_2D,[5 5],0,'both');
occ_2D_pad                        = flipud(occ_2D_pad);
occ_2D_pad_gau                    = zeros(33,33);


for idx_rc                        = 0:1:33*33-1
    idx_r                         = floor(idx_rc/33)+1;
    idx_c                         = mod(idx_rc,33)+1;
    temp                          = zeros(33,33);
    temp(idx_r,idx_c)             = occ_2D_pad(idx_r,idx_c);
    Sum_T                         = sum(temp(:));
    temp                          = imgaussfilt(temp, Sigma/2);
    temp                          = temp .* nbins_cir_2d_2d_pad;
    Sum_T2                        = sum(temp(:));
    Scale                         = Sum_T/Sum_T2;
    Scale(isnan(Scale))           = 1;
    temp                          = temp * Scale;
    Sum_T3                        = sum(temp(:));
    occ_2D_pad_gau                = occ_2D_pad_gau + temp;
end

occ_2D_pad_gau_depad              = occ_2D_pad_gau((1+5):(33-5), (1+5):(33-5));
occ_2D_pad_gau_depad(occ_2D_pad_gau_depad<occThresh)=0;
%----------Occupancy time
%--------------------------------------------------------------------------
 parfor idxnumShf                 = 1:numShuffles
        
    info                          = nan(ns);%% number of shuffles X number of cells   
     
    indexs                        = indexp+shift(idxnumShf,1);
    indexs2                       = mod(indexs, NoFrames);
    indexs2(indexs2==0)           = NoFrames;
    
    posx2                         = posx( indexs2 );
    posy2                         = posy( indexs2 );
            
    [~, ~, ~, xbin2, ybin2]       = histcounts2(posx2, posy2, bin_edges, bin_edges);
    % ---------------------------------------------------
    bin2_1D                       = (xbin2-1)*(size(bin_edges,2)-1)+ybin2;% Change the 2D bin indexes to 1D indexes for quick processing
    % ---------------------------------------------------    
    temp                          = zeros(length(occ), ns);%% Bin No. X cells number
   
    for segs                      = 1:ns
        spk_norm                  = spks(segs, Speed_r_th_idx)';
        temp2                     = accumarray( bin2_1D,spk_norm );%在不同bin里面每个cells激发次数？？？
        temp(1:length(temp2),segs)= temp2;
        disp([ 'Seg No:' num2str(segs) '/' num2str(ns) '--Shuffle No:' num2str(idxnumShf) '/' num2str(numShuffles)]);
    end
    
    temp_2D                       = reshape(temp,[(length(bin_edges)-1),(length(bin_edges)-1), ns]);
    temp_2D_pad                   = padarray(temp_2D,[5 5],0,'both');
    temp_2D_pad                   = flipud(temp_2D_pad);
    
    Spk_2D_pad_gau                = zeros(33,33,ns);
    
    for idx_rc                    = 0:1:33*33-1
        idx_r                     = floor(idx_rc/33)+1;
        idx_c                     = mod(idx_rc,33)+1;
        
        temp                      = zeros(33,33,ns);
        temp(idx_r,idx_c,:)       = temp_2D_pad(idx_r,idx_c,:);
        Sum_T                     = sum(temp,[1 2]);
        temp                      = imgaussfilt(temp, Sigma/2);
        temp                      = temp .* nbins_cir_2d_2d_pad;
        Sum_T2                    = sum(temp,[1 2]);
        Scale                     = Sum_T./Sum_T2;
        Scale(isnan(Scale))       = 1;
        temp                      = temp .* Scale;
        Sum_T3                    = sum(temp,[1 2]);
        Spk_2D_pad_gau            = Spk_2D_pad_gau + temp;
        
        disp([num2str(idx_r) 'X' num2str(idx_c) '/' num2str(33) 'X' num2str(33) '--Shuffle No:' num2str(idxnumShf) '/' num2str(numShuffles)]);
    end
    
    
    Spk_2D_pad_gau_depad          = Spk_2D_pad_gau((1+5):(33-5), (1+5):(33-5),:);
    
    FR_3D_gau                     = Spk_2D_pad_gau_depad./repmat(occ_2D_pad_gau_depad,1,1,ns);
    FR_3D_gau(FR_3D_gau==Inf)     = nan;
    
    FR_4D_gau(:,:,:,idxnumShf)    = FR_3D_gau;
    
    FR_3D_gau(isnan(FR_3D_gau))   = 0;    
    FR_3D_gau(FR_3D_gau==0)       = 0.000000000001;
    
    OCC_2D_gau                    = occ_2D_pad_gau_depad;
    
    OCC_2D_gau(isnan(OCC_2D_gau)) = 0;
    OCC_2D_gau(OCC_2D_gau==0)     = 0.000000000001;
    %---------------------------------
    P_occ                         = OCC_2D_gau./repmat(sum(OCC_2D_gau,[1 2]),size(OCC_2D_gau,1),size(OCC_2D_gau,2),1);
    refFR_mean                    = repmat(sum(repmat(P_occ,1,1,ns).*FR_3D_gau,[1 2]),size(OCC_2D_gau,1),size(OCC_2D_gau,2),1);
    refFR                         = FR_3D_gau./refFR_mean;
    %---------------------------------
    info                          = sum(P_occ.*refFR.*log2(refFR),[1 2]);%%log2(size(refFR,1)) is what?????
    info                          = abs(squeeze(info));
        
    info_Shuffle(idxnumShf,:)     = info;

 end

info_Real                         = info_Shuffle(1,:);
infoP                             = sum(repmat(info_Real,numShuffles-1,1) >= info_Shuffle(2:end,:),1)/(numShuffles-1);
display(['Percent over 0.95: ' num2str(100*sum(infoP>=0.95)/ns) '%']); %百分位在95%（稳定性在95%）的神经细胞占总细胞个数的百分比

end