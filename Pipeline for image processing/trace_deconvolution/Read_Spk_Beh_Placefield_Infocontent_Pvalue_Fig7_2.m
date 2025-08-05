%##########################################################################
clear all;
close all;
clc;
Folder             = 'D:\Cscope Data\PP367\Session1_12_51_06\Processed\';
addpath('D:\Cscope Data\PP367\Session1_12_51_06\Processed\Matlab Codes\Miniscope Analysis Software\Miniscope Analysis10-6-17\Miniscope Analysis\Util\','-begin');
%##########################################################################
binSize            = 2;        % 2 cm bins for an arena of 80cm diameter
binRad             = 23;
bin_edges          = [-binRad:binSize:binRad]; 
min_Entries        = 3; % must visit a bin at least 3 tims to be counted
occThresh          = 0.1; % 0.05s is 50ms
%##########################################################################
%Load deconvolved spikes from OASIS jupyter notebook
load([Folder 'Spk_Daniel_g.mat']);
spks               = s_2d;

figure(1);
for idx=1:20
plot(spks(idx,:)+80*idx)
hold on;
end
spks_binary        = spks>0;%Binarize the spikes bigger than zero.
%==========================================================================
Folder_Be = 'D:\Cscope Data\PP367\Session1_12_51_06\BehavCam_0_Analysis\Matlab Codes\Behav_Data_Prod';
load([Folder_Be '\behav_g1.mat']);
Start              = 1;%% start frame number for analysis

T_beh              = behav.time;
T_beh              = T_beh/1000;
T_beh              = T_beh(Start:end);%% timestamp from behaviral recording

dt                 = [0; diff(T_beh)];
dt_m               = median(dt);

posx               = behav.position(:,1);
posy               = behav.position(:,2);
posx               = posx(Start:end);
posy               = posy(Start:end); %% positions from behavioral recording

posx               = smoothts(posx','b',ceil(1/dt_m))';%Positions are smoothed around every second!
posy               = smoothts(posy','b',ceil(1/dt_m))';%Positions are smoothed around every second!

Speed              = behav.speed;
Speed              = Speed(Start:end); %% speed from behavioral recording
                  
dx                 = [0; diff(posx)];
dy                 = [0; diff(posy)];

Speed_r            = sqrt(dx.^2+dy.^2)./dt_m; %% use the median timestmap for speed!
Speed_r            = abs(smoothts(Speed_r','b',ceil(1/dt_m)))'; % need to use ' to make the matrix a column vector, not row vector

figure(100);
    plot(Speed,'r');
    hold on;
    plot(Speed_r,'k','LineWidth', 1.5)
title(['Median Speed is: ' num2str(median(Speed_r))])
%==========================================================================
T                  = readtable('D:\Cscope Data\PP367\Session1_12_51_06\Miniscope\timeStamps.csv','Range','B2:B43351');
T_spk              = table2array(T); %% timestamp from neural recording
T_spk              = T_spk/1000;
Dt_spk             = [0; diff(T_spk)];
%==========================================================================
T_Beh_Spk          = [T_beh; T_spk];
T_Beh_Spk          = sort(T_Beh_Spk);
% concatenate the timestamp from behavior and timestamp from neuron
%==========================================================================
posx               = interp1(T_beh,posx,T_Beh_Spk,'nearest');
posx               = posx-(min(posx)+max(posx))/2;
posx(posx>binRad)  = binRad;
posx(posx<-binRad) = -binRad;
posy               = interp1(T_beh,posy,T_Beh_Spk,'nearest');
posy               = posy-(min(posy)+max(posy))/2;
posy(posy>binRad)  = binRad;
posy(posy<-binRad) = -binRad;
[Flg,Spk_idx]      = ismember(T_spk,T_Beh_Spk,'rows');
posx               = posx(Spk_idx);
posy               = posy(Spk_idx);

figure(200)
    plot(posx)
    hold on;
    plot(posy)
%==========================================================================
End                = length(posx);
Start              = 1;
%==========================================================================
Speed_r            = interp1(T_beh,Speed_r,T_Beh_Spk,'nearest');
Speed_r            = Speed_r(Spk_idx);

Speed_r_th         = 5;  %Speed threshold is Speed_r_th cm/s
Speed_r_th_h       = 50; %Speed threshold is Speed_r_th cm/s
Speed_r_th_idx     = ((Speed_r>=Speed_r_th)&(Speed_r<=Speed_r_th_h))'; %speed > Speed_r_th cm/s
figure(300);
plot(Speed_r)            %Set the speed value to be true when it is between 5cm/s and 50cm/s
ylim([0 30])
figure(303);
h = histogram(Speed_r)
%==========================================================================
posx               = posx( Start:End );
posy               = posy( Start:End );
spks               = spks(:, Start:End );

spks_binary        = spks_binary(:, Start:End );
Speed_r            = Speed_r( Start:End );
Speed_r_th_idx     = Speed_r_th_idx( Start:End );
Dt_spk             = Dt_spk( Start:End );
Dt_spk             = Dt_spk';
T_spk              = T_spk( Start:End );
figure(400)
% plot(posx,posy,'r.')
xlim([-binRad binRad])
ylim([-binRad binRad])
xticks(-binRad:2:binRad)
yticks(-binRad:2:binRad)
grid on;
hold on;
plot(posx(Speed_r_th_idx'),posy(Speed_r_th_idx'),'b.');
%==========================================================================
ns                           = size(spks,1); %% number of cells
Sigma                        = 3;
%==========================================================================
    indexp                   = 1:1:length(T_spk);   
    indexs                   = indexp;
    indexs2                  = mod(indexs, length(T_spk));
    indexs2(indexs2==0)      = length(T_spk);  
    dt_m                     = median(Dt_spk);    
    occ                      = zeros((length(bin_edges)-1)*(length(bin_edges)-1),1); %% Bin No. X 1    
    [~, ~, ~, xbin, ybin]    = histcounts2(posx, posy, bin_edges, bin_edges);
    % ---------------------------------------------------
    bin2_1D                  = (xbin-1)*(size(bin_edges,2)-1)+ybin;% Change the 2D bin indexes to 1D indexes for quick processing
    % ---------------------------------------------------
    Dt_spk                   = Dt_spk';        
    
    CirX                     = -binRad:0.1:binRad;
    CirY                     = -binRad:0.1:binRad;
    [Cols Rows]              = meshgrid(CirX, CirY);
    centerX                  = 0;
    centerY                  = 0;
    Radius                   = binRad;
    Cols_1D                  = reshape(Cols,[],1);  
    Rows_1D                  = reshape(Rows,[],1);    
    Mask_Cir_2d              = (Cols - centerY).^1 + (Rows - centerX).^1 <= Radius.^1;
    Mask_Cir                 = (abs(Cols_1D - centerY) <= Radius)|(abs(Rows_1D - centerX)<= Radius);
    Cols_P                   = Cols(Mask_Cir);
    Rows_P                   = Rows(Mask_Cir);
    
    [~, ~, ~, xbin_cir, ybin_cir] = histcounts2(Cols_P, Rows_P, bin_edges, bin_edges);
    bin2_cir_1D              = (xbin_cir-1)*(size(bin_edges,2)-1)+ybin_cir;
    nbins_cir                = accumarray(bin2_cir_1D,1);
    nbins_cir_2d             = zeros((length(bin_edges)-1)*(length(bin_edges)-1),1); %% Bin No. X 1
    nbins_cir_2d(1:length(nbins_cir))= nbins_cir;
    nbins_cir_2d_2d          = reshape(nbins_cir_2d,[(length(bin_edges)-1),(length(bin_edges)-1)]);
    nbins_cir_2d_2d          = double(nbins_cir_2d_2d>0);
    figure(500)
    heatmap((nbins_cir_2d_2d));
    title('Circular mask for the shape of open field')
    nbins_cir_2d_2d_pad      = padarray(nbins_cir_2d_2d,[5 5],0,'both');
    figure(501)
    heatmap((nbins_cir_2d_2d_pad));
    title('Circular mask for the shape of open field after zeropadding')
%==========================================================================
%%-----------Calculating the entering times threshold for spikes to remove
bin2_1D                      = (xbin-1)*(size(bin_edges,2)-1)+ybin;%
bin2_1D_Test_diff            = [(diff(bin2_1D)); 0];
bin2_1D_Test_diff            = bin2_1D_Test_diff .*(Speed_r_th_idx');
bin2_1D_Test_diff_idx        = (bin2_1D_Test_diff~=0);
bin2_1D_thrd_minEntr         = bin2_1D.*bin2_1D_Test_diff_idx;
bin2_1D_thrd_minEntr         = bin2_1D_thrd_minEntr(bin2_1D_thrd_minEntr~=0);
bin2_1D_thrd_minEntr_T       = accumarray(bin2_1D_thrd_minEntr,1);
bin2_1D_thrd_minEntr_T_idx   = find(bin2_1D_thrd_minEntr_T>=min_Entries);
bin2_1D_MinEntr              = zeros(size(bin2_1D,1),1);

for idx=1:size(bin2_1D_thrd_minEntr_T_idx,1)
    bin2_1D_MinEntr(bin2_1D==bin2_1D_thrd_minEntr_T_idx(idx))=1;
end

figure(599)
plot(spks(8,:))
figure(598)

plot(posx,posy,'k','LineWidth' ,1);
xlim([-binRad binRad])
ylim([-binRad binRad])
xticks(-binRad:5:binRad)
yticks(-binRad:5:binRad)
axis equal
axis tight
% 
% idx=8
% v = VideoWriter([Folder 'PosTrace_Spk_Speed_time.avi'],'Motion JPEG AVI');
% v.Quality = 95;
% open(v)
% idx_spk=(spks_binary(8,:)&Speed_r_th_idx);
% for idx = 1:size(posx,1)
%     clear gca
%     gcf=figure(600)
%     
%     title([num2str(idx)]);
%     plot(posx(1:idx),posy(1:idx),'k','LineWidth' ,1);
%     axis equal
%     
%     hold on;
%     plot(posx(idx_spk(1:idx)),posy(idx_spk(1:idx)),'.r','MarkerSize',30)    
%     %      set(gca,'Visible','off')
%     set(gcf,'color','w');
%     
%     drawnow;
%   
%     set(gca, 'Xlim',[-binRad binRad])
%     set(gca, 'Ylim',[-binRad binRad])
%     xticks(-binRad:5:binRad)
%     yticks(-binRad:5:binRad)
%     axis equal
%     axis tight
%     grid on;
%     imagewd = getframe(gcf);
%     writeVideo(v,imagewd);
%     hold off
%     idx
% end
% close(v)

for idx=8
    gcf=figure(600)
        plot(posx,posy,'k','LineWidth' ,1);
        xlim([-binRad binRad])
        ylim([-binRad binRad])  
        xticks(-binRad:5:binRad)
        yticks(-binRad:5:binRad)
        axis equal
        axis tight
        grid on;
        hold on;
        plot(posx(spks_binary(idx,:)&Speed_r_th_idx),posy(spks_binary(idx,:)>0&Speed_r_th_idx),'.r','MarkerSize',30)
        title([num2str(idx)]);
    
        set(gca,'Visible','off')
        set(gcf,'color','w');
        hold off;
    imagewd = getframe(gcf); 
    imwrite(imagewd.cdata, [Folder 'PosTrace_Spk_Speed.tif'],'WriteMode','append');
    idx
end

% t= 1:1:size(posx,1);
% tt=[t;t];
% xx=[posx';posx'];
% yy=[posy';posy'];
% 
% figure(601)
% hs=surf(xx,yy,tt,tt,'EdgeColor','interp') %// color binded to "t" values
% clrmap = 'gray';
% title(clrmap)
% set(hs,'LineWidth',1)
% colormap(clrmap)
% view(2) %// view(0,90)
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
Htmap_flg                             = 1;%% only calculate the heatmap and info content

occ                                   = zeros((length(bin_edges)-1)*(length(bin_edges)-1),1); %% Bin No. X 1
FR_2D                                 = zeros(40,40,ns);

    info                              = nan(ns);%% number of shuffles X number of cells   
    [~, ~, ~, xbin, ybin]             = histcounts2(posx, posy, bin_edges, bin_edges);
    bin2_1D                           = (xbin-1)*(size(bin_edges,2)-1)+ybin;% Change the 2D bin indexes to 1D indexes for quick processing
  
Comment_flg                           = 0;    
if ~Comment_flg
    
    Dt_spk_spd                        = Dt_spk(Speed_r_th_idx');
    temp                              = accumarray(bin2_1D(Speed_r_th_idx),median(Dt_spk));
    occ(1:length(temp))               = temp; %% occupantion time in each bin!!!    
    occ_NoThrd                        = occ;
    occ_NoThrd(isnan(occ_NoThrd))                   = 0;
   
    occ(occ<occThresh)                = nan;    
    occ(isnan(occ))                   = 0;    
       
    occ_2D                            = reshape(occ',[(length(bin_edges)-1),(length(bin_edges)-1)]);
       
    occ_2D_Nothrd                     = reshape(occ_NoThrd',[(length(bin_edges)-1),(length(bin_edges)-1)]);
    occ_2D_pad                        = padarray(occ_2D_Nothrd,[5 5],0,'both');
    occ_2D_pad                        = flipud(occ_2D_pad);
    occ_2D_pad_gau                    = zeros(33,33);        
    
    for idx_rc                = 0:1:33*33-1
        idx_r                 = floor(idx_rc/33)+1;
        idx_c                 = mod(idx_rc,33)+1;
        temp                  = zeros(33,33);
        temp(idx_r,idx_c)     = occ_2D_pad(idx_r,idx_c);
         Sum_T                 = sum(temp(:));
        temp                  = imgaussfilt(temp, Sigma/2);
        temp                  = temp .* nbins_cir_2d_2d_pad;
        Sum_T2                = sum(temp(:));
        Scale                 = Sum_T/Sum_T2;
        Scale(isnan(Scale))   = 1;
        temp                  = temp * Scale;
        Sum_T3                = sum(temp(:));
        occ_2D_pad_gau        = occ_2D_pad_gau + temp;
    end    

    occ_2D_pad_nogau_depad    = flipud(occ_2D_Nothrd);
    occ_2D_pad_gau_depad      = occ_2D_pad_gau((1+5):(33-5), (1+5):(33-5));
    occ_2D_pad_gau_depad_nothrd = occ_2D_pad_gau_depad;
    occ_2D_pad_gau_depad(occ_2D_pad_gau_depad<occThresh)=0;
end    
    
    temp                           = zeros(length(occ), ns);%% Bin No. X cells number
    temp_b                         = zeros(length(occ), ns);%% Bin No. X cells number
    for segs                       = 1:ns
        spk_norm                   = spks(segs, Speed_r_th_idx)';
        temp2                      = accumarray( bin2_1D(Speed_r_th_idx),spk_norm );%在不同bin里面每个cells激发次数？？？
        temp(1:length(temp2),segs) = temp2;
        
        spk_norm                   = spks_binary(segs, Speed_r_th_idx)';
        temp2                      = accumarray( bin2_1D(Speed_r_th_idx),spk_norm );%在不同bin里面每个cells激发次数？？？
        temp_b(1:length(temp2),segs) = temp2;
        
        disp([ 'Seg No:' num2str(segs) '/' num2str(ns)]);
    end
    temp_2D                        = reshape(temp,[(length(bin_edges)-1),(length(bin_edges)-1), ns]);
    temp_2D_pad                    = padarray(temp_2D,[5 5],0,'both');
    temp_2D_pad                    = flipud(temp_2D_pad);        
    
    temp_b_2D                      = reshape(temp_b,[(length(bin_edges)-1),(length(bin_edges)-1), ns]);
    temp_b_2D_pad                  = padarray(temp_b_2D,[5 5],0,'both');
    temp_b_2D_pad                  = flipud(temp_b_2D_pad);  
   
    Spk_2D_pad_gau                 = zeros(33,33,ns);
    Spk_b_2D_pad_gau               = zeros(33,33,ns);
    
    for idx_rc             = 0:1:33*33-1
        idx_r                 = floor(idx_rc/33)+1;
        idx_c                 = mod(idx_rc,33)+1;
       
        temp                   = zeros(33,33,ns);
        temp(idx_r,idx_c,:)    = temp_2D_pad(idx_r,idx_c,:);            
        Sum_T                  = sum(temp,[1 2]);
        temp                   = imgaussfilt(temp, Sigma/2);
        temp                   = temp .* nbins_cir_2d_2d_pad;
        Sum_T2                 = sum(temp,[1 2]);
        Scale                  = Sum_T./Sum_T2;
        Scale(isnan(Scale))    = 1;
        temp                   = temp .* Scale;
        Sum_T3                 = sum(temp,[1 2]);
        Spk_2D_pad_gau         = Spk_2D_pad_gau + temp;
        
        temp                   = zeros(33,33,ns);
        temp(idx_r,idx_c,:)    = temp_b_2D_pad(idx_r,idx_c,:);            
        Sum_T                  = sum(temp,[1 2]);
        temp                   = imgaussfilt(temp, Sigma/2);
        temp                   = temp .* nbins_cir_2d_2d_pad;
        Sum_T2                 = sum(temp,[1 2]);
        Scale                  = Sum_T./Sum_T2;
        Scale(isnan(Scale))    = 1;
        temp                   = temp .* Scale;
        Sum_T3                 = sum(temp,[1 2]);
        Spk_b_2D_pad_gau       = Spk_b_2D_pad_gau + temp;
        
        disp([num2str(idx_r) 'X' num2str(idx_c) '/' num2str(45) 'X' num2str(45)]);
    end
    
    Spk_2D_pad_nogau_depad     = flipud(temp_2D);
    Spk_b_2D_pad_nogau_depad   = flipud(temp_b_2D);
    Spk_2D_pad_gau_depad       = Spk_2D_pad_gau((1+5):(33-5), (1+5):(33-5),:);
    Spk_b_2D_pad_gau_depad     = Spk_b_2D_pad_gau((1+5):(33-5), (1+5):(33-5),:);   
    
    FR_3D_gau                          = Spk_2D_pad_gau_depad./repmat(occ_2D_pad_gau_depad,1,1,ns);
    FR_3D_gau(FR_3D_gau==Inf)          = nan;
    FR_3D_gau(isnan(FR_3D_gau))        = 0;

    FR_3D_gau(FR_3D_gau==0)            = 0.000000000001;
    OCC_2D_gau = occ_2D_pad_gau_depad;
    
    OCC_2D_gau(isnan(OCC_2D_gau))              = 0;
    OCC_2D_gau(OCC_2D_gau==0)          = 0.000000000001;    
    
    %---------------------------------
    P_occ                              = OCC_2D_gau./repmat(sum(OCC_2D_gau,[1 2]),size(OCC_2D_gau,1),size(OCC_2D_gau,2),1);
    refFR_mean                         = repmat(sum(repmat(P_occ,1,1,ns).*FR_3D_gau,[1 2]),size(OCC_2D_gau,1),size(OCC_2D_gau,2),1);
    refFR                              = FR_3D_gau./refFR_mean;
    %---------------------------------
    info                               = sum(P_occ.*refFR.*log2(refFR),[1 2]);%%log2(size(refFR,1)) is what?????
    info                               = abs(squeeze(info));       
       
    for segs                           = 1:ns
        
        gcf=figure(1001)
        gcf.Position                   = [100 100 800 600];
        ClrBar_Font                    = 10;
        ha = tight_subplot(3,3,[.04 .0],[.02 .1],[.02 .02])
        sgtitle(['Cell: #' num2str(segs) '(occu thresd:' num2str(occThresh) 's)']);
        
        Spk_2D_pad_nogau_temp          = Spk_2D_pad_nogau_depad(:,:,segs);
        Spk_2D_pad_gau_depad_temp      = Spk_2D_pad_gau_depad(:,:,segs);
        Spk_b_2D_pad_nogau_depad_temp  = Spk_b_2D_pad_nogau_depad(:,:,segs);
        Spk_b_2D_pad_gau_depad_temp    = Spk_b_2D_pad_gau_depad(:,:,segs);
        
        axes(ha(3))
        cm3                            = colormap(ha(3),'jet');
        pos_3                          = get(ha(3),'position');
        pos_3(1)                       = pos_3(1);
        subplot('Position',pos_3);
        plot(posx,posy,'k','LineWidth' ,0.5);
        xlim([-binRad binRad])
        ylim([-binRad binRad])
        xticks(-binRad:5:binRad)
        yticks(-binRad:5:binRad)
        axis equal
        axis tight
        hold on;
        posx_spk_spd                   = posx(spks_binary(segs,:)&Speed_r_th_idx);
        posy_spk_spd                   = posy(spks_binary(segs,:)&Speed_r_th_idx);
        spks_spk_spd                   = spks(segs,spks_binary(segs,:)&Speed_r_th_idx)';
        spk_color                      = floor(spks_spk_spd/max(spks_spk_spd(:))*255)+1;
        scatter(posx_spk_spd,posy_spk_spd, 20,cm3(spk_color,:) ,'filled');
        colormap('jet');
        c3                             = colorbar;
        title(c3,'P','FontSize',ClrBar_Font);
        axis off;
        set(gcf,'color','w');
        
        axes(ha(1))
        imagesc(occ_2D_pad_nogau_depad);
        caxis([0 5]);
        cm1                            = colormap(ha(1), 'parula')
        c1                             = colorbar;
        title(c1,'S','FontSize',ClrBar_Font);
        title('Occu time(s):No Gau');
        axis equal;
        xlim([1 40])
        ylim([1 40])
        axis tight
        axis off
        
        axes(ha(2))
        imagesc(occ_2D_pad_gau_depad_nothrd);
        caxis([0 3]);
        cm2                           = colormap(ha(2), 'parula')
        c2                            = colorbar;
        title(c2,'S','FontSize',ClrBar_Font);
        title(['Occu time(s):\sigma=' num2str(Sigma)]);
        axis equal;
        xlim([1 40])
        ylim([1 40])
        axis tight
        axis off
        
        axes(ha(4))
        imagesc(Spk_2D_pad_nogau_temp);
        Max_V                        = max(Spk_2D_pad_nogau_temp(:));
        Max_V                        = (floor(Max_V/2)+1)*2;
        caxis([0 Max_V]);
        cm4                          = colormap(ha(4), 'parula')
        c4                           = colorbar;
        title(c4,'V','FontSize',ClrBar_Font);
        title( ['Neu Act:No Gau']);
        axis equal;
        xlim([1 40])
        ylim([1 40])
        axis tight
        axis off
        
        axes(ha(5))
        imagesc(Spk_2D_pad_gau_depad_temp);
        Max_V                       = max(Spk_2D_pad_gau_depad_temp(:));
        Max_V                       = (floor(Max_V/2)+1)*2;
        caxis([0 Max_V]);
        cm5                         = colormap(ha(5), 'parula')
        axis equal;
        xlim([1 40])
        ylim([1 40])
        axis tight
        axis off
        c5                          = colorbar;
        title(c5,'V','FontSize',ClrBar_Font);
        title(['Neu Act:\sigma=' num2str(Sigma)]);
        
        axes(ha(7))
        imagesc(Spk_b_2D_pad_nogau_depad_temp);
        Max_V                       = max(Spk_b_2D_pad_nogau_depad_temp(:));
        Max_V                       = (floor(Max_V/2)+1)*2;
        caxis([0 Max_V]);
        cm7                         = colormap(ha(7), 'parula')
        axis equal;
        xlim([1 40])
        ylim([1 40])
        axis tight
        axis off
        c7                          = colorbar;
        title(c7,'F','FontSize',ClrBar_Font);
        title( ['Neu Act Binary:No Gau']);
        
        axes(ha(8))
        imagesc(Spk_b_2D_pad_gau_depad_temp);
        Max_V                       = max(Spk_b_2D_pad_gau_depad_temp(:));
        Max_V                       = (floor(Max_V/2)+1)*2;
        caxis([0 Max_V]);
        cm8                         = colormap(ha(8), 'parula')
        axis equal;
        xlim([1 40])
        ylim([1 40])
        axis tight
        axis off
        c8                          = colorbar;
        title(c8,'F','FontSize',ClrBar_Font);
        title(['Neu Act Binary:\sigma=' num2str(Sigma)]);
        
        axes(ha(6))
        FR_2D_gau                   = Spk_2D_pad_gau_depad_temp./occ_2D_pad_gau_depad;
        FR_2D_gau(FR_2D_gau==Inf)   = nan;
        imAlpha=ones(size(FR_2D_gau));
        imAlpha(isnan(FR_2D_gau))=0;
        imagesc(FR_2D_gau,'AlphaData',imAlpha);
        Max_V                       = max(FR_2D_gau(:));
        Max_V                       = (floor(Max_V/2)+1)*2;
        caxis([0 Max_V]);
        cm6                         = colormap(ha(6), 'parula')
        axis equal;
        xlim([1 40])
        ylim([1 40])
        axis tight
        axis off
        c6                          = colorbar;
        title(c6,'R','FontSize',ClrBar_Font);
        title(['FR Htmp: Info=' num2str(info(segs))]);%
        
        axes(ha(9))
        FR_b_2D_gau                   = Spk_b_2D_pad_gau_depad_temp./occ_2D_pad_gau_depad;
        FR_b_2D_gau(FR_b_2D_gau==Inf) = nan;
        imAlpha=ones(size(FR_b_2D_gau));
        imAlpha(isnan(FR_b_2D_gau))=0;
        imagesc(FR_b_2D_gau,'AlphaData',imAlpha);
        Max_V                         = max(FR_b_2D_gau(:));
        Max_V                         = (floor(Max_V/2)+1)*2;
        caxis([0 Max_V]);
        cm9                           = colormap(ha(9), 'parula')
        axis equal;
        xlim([1 40])
        ylim([1 40])
        axis tight
        axis off
        c9                            = colorbar;
        title(c9,'R','FontSize',ClrBar_Font);
        title(['FR Htmp Binary:\sigma=' num2str(Sigma)]);
        
        drawnow;
        
        imagewd                       = getframe(gcf);
        imwrite(imagewd.cdata, [Folder 'Occ_Neural_' num2str(Sigma) '_2.tif'],'WriteMode','append');
        disp([num2str(segs),'/' num2str(ns)]);
        
        close(gcf)
    end
% %--------------------------------------------------------------------------
numShuffles                       = 500;
% [Spk_2D_pad_gau_depad_80,occ_2D_pad_gau_depad, info_Shuffle2, infoP, shift]        = Info_Gau_ClG_Fig7_2_cell(posx, posy, Speed_r_th_idx, Dt_spk,occThresh, bin_edges, spks, Sigma, numShuffles);
% [info_Shuffle,infoP] = Info_Gau_ClG_Fig7_2(posx, posy, Speed_r_th_idx, Dt_spk,occThresh, bin_edges, spks, Sigma, numShuffles);
[FR_4D_gau, info_Shuffle,infoP] = Info_Gau_ClG_Shuffle1st_Placefield_Fig7(posx, posy, Speed_r_th_idx, Dt_spk,occThresh, bin_edges, spks, Sigma, numShuffles);
% %--------------------------------------------------------------------------
save([Folder 'Info_P_FRshlf.mat'],'info_Shuffle','infoP','FR_4D_gau','-v7.3');
%%--------------------------------------------------------------------------    
load([Folder 'Info_P_FRshlf.mat']);
size(FR_4D_gau)
Heatmaps_seg_1D_Shlf=[];
for segs=1:ns
    
    Heatmaps_seg = squeeze(FR_4D_gau(:,:,segs,:));
    Heatmaps_seg_1D = reshape(Heatmaps_seg, [1 , size(Heatmaps_seg,1)*size(Heatmaps_seg,2)*size(Heatmaps_seg,3)]);
    Heatmaps_1D(segs,:) = Heatmaps_seg_1D(~isnan(Heatmaps_seg_1D));
    Heatmaps95(segs) = prctile(Heatmaps_1D(segs,:),95);
    segs
end


figure(9003);
plot(Heatmaps_1D(3,:))
figure(9004);
histogram(Heatmaps_1D(3,:))
sum(Heatmaps_1D(1,:)>=51.599933928571350)
sum(Heatmaps_1D(1,:)<51.599933928571350)/615228

Neuron_4     = [3];
figure(9005)
for idx=1:length(Neuron_4)
    
    Spk_2D_pad_gau_depad_temp      = Spk_2D_pad_gau_depad(:,:,Neuron_4(idx));
    FR_2D_gau                      = Spk_2D_pad_gau_depad_temp./occ_2D_pad_gau_depad;
    FR_2D_gau(FR_2D_gau==Inf)   = nan;
    imAlpha=ones(size(FR_2D_gau));
    imAlpha(isnan(FR_2D_gau))=0;
    imagesc(FR_2D_gau,'AlphaData',imAlpha);
    Max_V                       = max(FR_2D_gau(:));
    Max_V                       = (floor(Max_V/2)+1)*2;
    caxis([0 Max_V]);
    axis equal;
    xlim([1 40])
    ylim([1 40])
    axis tight
    axis off
    set(gcf,'color','w');
    %     c6                          = colorbar;
    %     title(c6,'R','FontSize',ClrBar_Font);
    %     title(['FR Htmp: Info=' num2str(info(segs))]);%
    hold on;
    FR_2D_gau                      = Spk_2D_pad_gau_depad_temp./occ_2D_pad_gau_depad;
    FR_2D_gau(FR_2D_gau==Inf)      = 0;
    FR_2D_gau_norm                 = FR_2D_gau/max(FR_2D_gau(:));
    FR_2D_gau_Binarized           = (FR_2D_gau_norm>=0.05);
    [B,L] = bwboundaries(FR_2D_gau_Binarized);           
%     for k=1:length(B)
%         b = B{k};
%         Area(k) = polyarea(b(:,2),b(:,1));
%         plot(b(:,2),b(:,1),'r','LineWidth',2);
%         hold on;
%     end
    hold on;
    
    FR_2D_gau_Bina_2              = (FR_2D_gau>=Heatmaps95(Neuron_4(idx)));    
    [B2,L] = bwboundaries(FR_2D_gau_Bina_2);      
%     for k=1:length(B2)
%         b2 = B2{k};
%         Area2(k) = polyarea(b2(:,2),b2(:,1));
%         plot(b2(:,2),b2(:,1),'g','LineWidth',2);
%         hold on;
%     end
    hold on;
    
    for B2idx = 1:length(B2)
         b2 = B2{B2idx};
        for B1idx = 1:length(B)
             b = B{B1idx};
    [in,on] = inpolygon(b2(1,2),b2(1,1),b(:,2),b(:,1));
     if((in)||(on))
        
        plot(b(:,2),b(:,1),'r','LineWidth',2);
        hold on;
    end
    hold on; 
         
         
     end
        
        end
        
       
   
end

Neurongood = [3,7,21,24,38,39,40,45,73,90,109,111,114,115,116,126,139,143,157,174,192,195,199,213,225];

ClorNo        = length(Neurongood);
clr           = lines( ClorNo );
Neuron_4     = [3,7,21,109];
[~,idx_cell] = ismember( Neuron_4 , Neurongood);
gcf=figure(9006)
gcf.Position                   = [100 100 800 800];
ClrBar_Font                    = 10;
ha = tight_subplot(4,4,[.02 .0],[.02 .1],[.02 .02])
sgtitle(['Threshold: ' num2str(0.5)]);
for idx=1:length(Neuron_4)
    
    axes(ha(idx))
    plot(posx,posy,'k','LineWidth' ,2);
    xlim([-40 40])
    ylim([-40 40])
    xticks(-40:5:40)
    yticks(-40:5:40)
    axis equal
    axis tight
    hold on;
    
    plot(posx(spks(Neuron_4(idx),:)>0),posy(spks(Neuron_4(idx),:)>0),'.','MarkerSize',20,'Color', clr(idx_cell(idx),:))
    set(gca,'Visible','off')
    set(gcf,'color','w');
    hold on;
    axes(ha(idx+4))
    Spk_2D_pad_gau_depad_temp      = Spk_2D_pad_gau_depad(:,:,Neuron_4(idx));
    FR_2D_gau                      = Spk_2D_pad_gau_depad_temp./occ_2D_pad_gau_depad;
    FR_2D_gau(FR_2D_gau==Inf)   = nan;
    imAlpha=ones(size(FR_2D_gau));
    imAlpha(isnan(FR_2D_gau))=0;
    imagesc(FR_2D_gau,'AlphaData',imAlpha);
    Max_V                       = max(FR_2D_gau(:));
    Max_V                       = (floor(Max_V/2)+1)*2;
    caxis([0 Max_V]);
    axis equal;
    xlim([1 40])
    ylim([1 40])
    axis tight
    axis off
    set(gcf,'color','w');
    %     c6                          = colorbar;
    %     title(c6,'R','FontSize',ClrBar_Font);
    %     title(['FR Htmp: Info=' num2str(info(segs))]);%
    hold on;
    FR_2D_gau_norm                 = FR_2D_gau/max(FR_2D_gau(:));
    FR_2D_gau_Binarized           = (FR_2D_gau_norm>=0.05);
    [B,L] = bwboundaries(FR_2D_gau_Binarized);
    for k=1:length(B)
        b = B{k};
        plot(b(:,2),b(:,1),'r','LineWidth',2);
        hold on;
    end
    hold on;
    
    axes(ha(idx+8))
    
    imagesc(FR_2D_gau,'AlphaData',imAlpha);
    Max_V                       = max(FR_2D_gau(:));
    Max_V                       = (floor(Max_V/2)+1)*2;
    caxis([0 Max_V]);
    axis equal;
    xlim([1 40])
    ylim([1 40])
    axis tight
    axis off
    set(gcf,'color','w');
    hold on;
    FR_2D_gau_Bina_2              = (FR_2D_gau>=Heatmaps95(Neuron_4(idx)));
    [B2,L] = bwboundaries(FR_2D_gau_Bina_2);
    FR_2D_gau_BW=zeros(size(FR_2D_gau,1),size(FR_2D_gau,1));
    for B2idx = 1:length(B2)
        b2 = B2{B2idx};
        for B1idx = 1:length(B)
            b = B{B1idx};
            [in,on] = inpolygon(b2(1,2),b2(1,1),b(:,2),b(:,1));
            if((in)||(on))
                plot(b(:,2),b(:,1),'r','LineWidth',2);
                hold on;
                
                
                Mask_temp = roipoly(FR_2D_gau,b(:,2),b(:,1));
                %                Mask_temp = poly2mask(b(:,2),b(:,1),size(FR_2D_gau,1),size(FR_2D_gau,1));
                for poidx = 1:size(b,1)
                    Mask_temp(b(poidx,1),b(poidx,2))=1;
                end
                FR_2D_gau_BW = FR_2D_gau_BW+Mask_temp;
                
            end
            hold on;
        end
    end
    axes(ha(idx+12))
    FR_2D_gau_BW = FR_2D_gau_BW>0;
    imagesc(FR_2D_gau_BW,'AlphaData',imAlpha);
    Max_V                       = max(FR_2D_gau_BW(:));
    Max_V                       = (floor(Max_V/2)+1)*2;
    caxis([0 Max_V]);
    axis equal;
    xlim([1 40])
    ylim([1 40])
    axis tight
    axis off
    set(gcf,'color','w');
    
    
end
drawnow;
hold off
%%-------------------------------------------------------------------------
%--------------------------------------------------------------------------
%%-------------------------------------------------------------------------
vw = VideoWriter(['D:\Cscope Data\PP367\Session1_12_51_06\Processed\' 'PlaceField_Check2.avi'],'Motion JPEG AVI');
vw.Quality = 95;
open(vw)
infoP_95 = find(infoP>0.95);

Placefield_Sum = zeros(size(FR_2D_gau,1),size(FR_2D_gau,1));

for Neuron_idx                        = 1:4:size(infoP_95,2)-4
    
    Neuron_4                          = infoP_95(Neuron_idx:Neuron_idx+3);
    gcf                               = figure(9006)
    gcf.Position                      = [100 100 800 600];
    ClrBar_Font                       = 10;
    ha                                = tight_subplot(3,4,[.02 .0],[.02 .1],[.02 .02])
    sgtitle(['Threshold: ' num2str(0.5)]);
    for idx                           = 1:length(Neuron_4)
                
        Spk_2D_pad_gau_depad_temp     = Spk_2D_pad_gau_depad(:,:,Neuron_4(idx));
        FR_2D_gau                     = Spk_2D_pad_gau_depad_temp./occ_2D_pad_gau_depad;
        FR_2D_gau(FR_2D_gau==Inf)     = nan;
        imAlpha                       = ones(size(FR_2D_gau));
        imAlpha(isnan(FR_2D_gau))     = 0;
        FR_2D_gau_norm                = FR_2D_gau/max(FR_2D_gau(:));
        FR_2D_gau_Binarized           = (FR_2D_gau_norm>=0.05);          
        
        axes(ha(idx))
        plot(posx,posy,'k','LineWidth' ,2);
        xlim([-40 40])
        ylim([-40 40])
        xticks(-40:5:40)
        yticks(-40:5:40)
        axis equal
        axis tight
        hold on;
        
        plot(posx(spks(Neuron_4(idx),:)>0),posy(spks(Neuron_4(idx),:)>0),'.','MarkerSize',20,'Color', 'r')
        set(gca,'Visible','off')
        set(gcf,'color','w');
        hold on;
        
        axes(ha(idx+4))        
        imagesc(FR_2D_gau,'AlphaData',imAlpha);
        Max_V                       = max(FR_2D_gau(:));
        Max_V                       = (floor(Max_V/2)+1)*2;
        caxis([0 Max_V]);
        axis equal;
        xlim([1 40])
        ylim([1 40])
        axis tight
        axis off
        set(gcf,'color','w');
        hold on;
        
        [B,L]                         = bwboundaries(FR_2D_gau_Binarized);
        FR_2D_gau_Bina_2              = (FR_2D_gau>=Heatmaps95(Neuron_4(idx)));
        [B2,L]                        = bwboundaries(FR_2D_gau_Bina_2);
        FR_2D_gau_BW                  = zeros(size(FR_2D_gau,1),size(FR_2D_gau,1));
        FR_2D_gau_BW2                 = zeros(size(FR_2D_gau,1),size(FR_2D_gau,1));      
        for Bidx = 1:length(B)
            b = B{Bidx};
            b2_size = 0;
            for B21idx = 1:length(B2)
                b2 = B2{B21idx};
                [in,on] = inpolygon(b2(1,2),b2(1,1),b(:,2),b(:,1));
                if((in)||(on))%((in)||(on))
                    plot(b(:,2),b(:,1),'r','LineWidth',0.5);
                    hold on;
                    plot(b2(:,2),b2(:,1),'g','LineWidth',1);
                    hold on;
                    if(size(b2,1)<=5)
                        plot(b2(:,2),b2(:,1),'.g','MarkerSize',10);
                    end
                    
                    if(size(b2,1)>5)
                        
                        Mask_temp  = roipoly(FR_2D_gau,b(:,2),b(:,1));
                        for poidx  = 1:size(b,1)
                            Mask_temp(b(poidx,1),b(poidx,2))=1;
                        end
                        Mask_temp2 = roipoly(FR_2D_gau,b2(:,2),b2(:,1));
                        for poidx  = 1:size(b2,1)
                            Mask_temp2(b2(poidx,1),b2(poidx,2))=1;
                        end
                        
                        FR_2D_gau_BW2 = FR_2D_gau_BW2 + Mask_temp2;
                        FR_2D_gau_BW  = FR_2D_gau_BW + Mask_temp;
                    end
                end
            end          
            
        end
        axes(ha(idx+8))
        FR_2D_gau_BW  = FR_2D_gau_BW>0;
        FR_2D_gau_BW2 = FR_2D_gau_BW2>0;
        FR_2D_gau_BW_12 = FR_2D_gau_BW2*3 +FR_2D_gau_BW;
        imagesc(FR_2D_gau_BW_12,'AlphaData',imAlpha);
        Max_V                       = max(FR_2D_gau_BW_12(:));
        Max_V                       = (floor(Max_V/2)+1)*2;
        caxis([0 Max_V]);
        axis equal;
        xlim([1 40])
        ylim([1 40])
        axis tight
        axis off
        set(gcf,'color','w');
        
        Placefield_Sum = Placefield_Sum + FR_2D_gau_BW;
        
    end
    drawnow;
    imagewd = getframe(gcf);
    writeVideo(vw,imagewd);
    hold off
    close(gcf)
    
    
end
close(vw)

gcf=figure(9007);
imagesc(Placefield_Sum,'AlphaData',imAlpha);
Max_V                       = max(Placefield_Sum(:));
Max_V                       = (floor(Max_V/2)+1)*2.5;
caxis([0 Max_V]);
axis equal;
xlim([1 40])
ylim([1 40])
axis tight
axis off
set(gcf,'color','w');

gcf=figure(9008);
Placefield_Sum_Per = Placefield_Sum/sum(infoP>=0.95)*100
imagesc(Placefield_Sum_Per,'AlphaData',imAlpha);
Max_V   =  max(Placefield_Sum_Per(:))*1.3;
c1=colorbar;
title(c1,'%','FontSize',13);
caxis([0 Max_V]);
axis equal;
xlim([1 40])
ylim([1 40])
axis tight
axis off
set(gcf,'color','w');
imagewd = getframe(gcf);
imwrite(imagewd.cdata, [Folder 'Sum_PlaceField_Figure7h.tif']);

%%-------------------------------------------------------------------------
%%-------------------------------------------------------------------------
%%-------------------------------------------------------------------------
%%-------------------------------------------------------------------------
%%-------------Plot traces and heatmaps of example cells
 Neurongood = [3,21,24,38,39,45,7,73,90,111,114,115,116,126,139,40,143,157,174,192,195,199,213,109,225];
    
 ClorNo        = length(Neurongood);
 clr           = lines( ClorNo );   
    
    
Dt_spk_m   = median(Dt_spk);

% Time_plot =0:Dt_spk_m:Dt_spk_m*(size(Dt_spk,1)-1);

NumFrame   = 7248;%size(Dt_spk,1);
NumFrame   = size(Dt_spk,1);
Time_plot =0:Dt_spk_m:Dt_spk_m*(NumFrame-1);

% Ctraces    = c_2d(:,Start:End);
Ctraces    = C_Raw_g_flt(:,Start:End);

for index=1:length(Neurongood)   
    gcf9=figure(9009)
    Ctraces_norm = Ctraces(Neurongood(index),:)/max(Ctraces(Neurongood(index),:));    
    plot(Time_plot,Ctraces_norm(1,1:NumFrame)+index*1,'Color', clr(index,:),'LineWidth' ,1.5)
%     set(gca,'Visible','off')
    set(gcf9,'color','w');
    box off;
    hold on;    
end  
imagewd = getframe(gcf9);
imwrite(imagewd.cdata, [Folder 'Traces_Chosen_Axisoff_Fig7i_' num2str(floor(Dt_spk_m*(NumFrame-1))) 's.tif']);
% exportgraphics(gcf9,[Folder 'Traces_Chosen_Axisoff_Fig7i_' num2str(floor(Dt_spk_m*(NumFrame-1))) 's.tif'],'Resolution',1200)


Neuron_4     = [3,7,40,109];
[~,idx_cell] = ismember( Neuron_4 , Neurongood); 
for idx=1:length(Neuron_4)
    
    gcf=figure(3000)
    plot(posx,posy,'k','LineWidth' ,2);
    xlim([-40 40])
    ylim([-40 40])  
    xticks(-40:5:40)
    yticks(-40:5:40)
    axis equal
    axis tight
    hold on;
    plot(posx(spks(Neuron_4(idx),:)>0),posy(spks(Neuron_4(idx),:)>0),'.','MarkerSize',20,'Color', clr(idx_cell(idx),:))
    title([num2str(Neuron_4(idx))]);
    
    set(gca,'Visible','off')
    set(gcf,'color','w');
    hold off;
    imagewd = getframe(gcf);
    imwrite(imagewd.cdata, [Folder 'Cells_SFR_Traces_Mark_Fig7g.tif'],'WriteMode','append');
    idx    

    Spk_2D_pad_gau_depad_temp     = Spk_2D_pad_gau_depad(:,:,Neuron_4(idx));
    FR_2D_gau                     = Spk_2D_pad_gau_depad_temp./occ_2D_pad_gau_depad;
    FR_2D_gau(FR_2D_gau==Inf)     = nan;
    imAlpha                       = ones(size(FR_2D_gau));
    imAlpha(isnan(FR_2D_gau))     = 0;
    FR_2D_gau_norm                = FR_2D_gau/max(FR_2D_gau(:));
    FR_2D_gau_Binarized           = (FR_2D_gau_norm>=0.05);    
    
    [B,L]                         = bwboundaries(FR_2D_gau_Binarized);
    FR_2D_gau_Bina_2              = (FR_2D_gau>=Heatmaps95(Neuron_4(idx)));
    [B2,L]                        = bwboundaries(FR_2D_gau_Bina_2);
    FR_2D_gau_BW                  = zeros(size(FR_2D_gau,1),size(FR_2D_gau,1));
    FR_2D_gau_BW2                 = zeros(size(FR_2D_gau,1),size(FR_2D_gau,1));
    
    
    gcf=figure(3001)
    imagesc(FR_2D_gau,'AlphaData',imAlpha);
    Max_V                       = max(FR_2D_gau(:));
    Max_V                       = (floor(Max_V/2)+1)*2;
    caxis([0 Max_V]);
    axis equal;
    xlim([1 40])
    ylim([1 40])
    axis tight
    axis off
    set(gcf,'color','w');
    hold on;
    for Bidx = 1:length(B)
        b = B{Bidx};
        b2_size = 0;
        for B21idx = 1:length(B2)
            b2 = B2{B21idx};
            [in,on] = inpolygon(b2(1,2),b2(1,1),b(:,2),b(:,1));
            if((in)||(on))
                if(size(b2,1)>5)
                    plot(b(:,2),b(:,1),'r','LineWidth',2);
                    hold on;
                end
            end
        end
    end
    
    
    drawnow;
        
    imagewd = getframe(gcf);
    imwrite(imagewd.cdata, [Folder 'Cells_Heatmaps_Fig7g.tif'],'WriteMode','append');   
    
end



gcf=figure(3001)
gcf.Position                   = [100 100 800 600];
ClrBar_Font                    = 10;
ha = tight_subplot(3,4,[.02 .0],[.02 .1],[.02 .02])
sgtitle(['Threshold: ' num2str(0.5)]);
for idx=1:length(Neuron_4)
     
    axes(ha(idx))
    plot(posx,posy,'k','LineWidth' ,2);
    xlim([-40 40])
    ylim([-40 40])  
    xticks(-40:5:40)
    yticks(-40:5:40)
    axis equal
    axis tight
    hold on;    
  
    plot(posx(spks(Neuron_4(idx),:)>0),posy(spks(Neuron_4(idx),:)>0),'.','MarkerSize',20,'Color', clr(idx_cell(idx),:))
    set(gca,'Visible','off')
    set(gcf,'color','w');
    hold on;
    axes(ha(idx+4))   
    Spk_2D_pad_gau_depad_temp      = Spk_2D_pad_gau_depad(:,:,Neuron_4(idx));
    FR_2D_gau                      = Spk_2D_pad_gau_depad_temp./occ_2D_pad_gau_depad;
    FR_2D_gau(FR_2D_gau==Inf)   = nan;
    imAlpha=ones(size(FR_2D_gau));
    imAlpha(isnan(FR_2D_gau))=0;
    imagesc(FR_2D_gau,'AlphaData',imAlpha);
    Max_V                       = max(FR_2D_gau(:));
    Max_V                       = (floor(Max_V/2)+1)*2;
    caxis([0 Max_V]);  
    axis equal;
    xlim([1 40])
    ylim([1 40])
    axis tight
    axis off
    set(gcf,'color','w');
    %     c6                          = colorbar;
    %     title(c6,'R','FontSize',ClrBar_Font);
    %     title(['FR Htmp: Info=' num2str(info(segs))]);%
    hold on;
    axes(ha(idx+8))  
    FR_2D_gau                      = Spk_2D_pad_gau_depad_temp./occ_2D_pad_gau_depad;
    FR_2D_gau(FR_2D_gau==Inf)      = 0;
    FR_2D_gau_norm                 = FR_2D_gau/max(FR_2D_gau(:));
    FR_2D_gau_Binarized           = (FR_2D_gau_norm>=0.5);
    imagesc(FR_2D_gau_Binarized,'AlphaData',imAlpha);
    Max_V                       = max(FR_2D_gau_Binarized(:));
    Max_V                       = (floor(Max_V/2)+1)*2;
    caxis([0 Max_V]);  
    axis equal;
    xlim([1 40])
    ylim([1 40])
    axis tight
    axis off
    set(gcf,'color','w');    
    hold on;          
        
end
drawnow;
imagewd = getframe(gcf);
imwrite(imagewd.cdata, [Folder 'Cells_Pos_Heatmap_Placefield.tif'],'WriteMode','append');



gcf=figure(3002)
gcf.Position                   = [100 100 800 200];
ClrBar_Font                    = 10;
ha = tight_subplot(2,4,[.02 .0],[.02 .1],[.02 .02])
sgtitle(['Threshold: ' num2str(0.5)]);
for idx=1:length(Neuron_4)
     
    axes(ha(idx))
    plot(posx,posy,'k','LineWidth' ,2);
    xlim([-40 40])
    ylim([-40 40])  
    xticks(-40:5:40)
    yticks(-40:5:40)
    axis equal
    axis tight
    hold on;    
  
    plot(posx(spks(Neuron_4(idx),:)>0),posy(spks(Neuron_4(idx),:)>0),'.','MarkerSize',20,'Color', clr(idx_cell(idx),:))
    set(gca,'Visible','off')
    set(gcf,'color','w');
    hold on;
    axes(ha(idx+4))   
    Spk_2D_pad_gau_depad_temp      = Spk_2D_pad_gau_depad(:,:,Neuron_4(idx));
    FR_2D_gau                      = Spk_2D_pad_gau_depad_temp./occ_2D_pad_gau_depad;
    FR_2D_gau(FR_2D_gau==Inf)   = nan;
    imAlpha=ones(size(FR_2D_gau));
    imAlpha(isnan(FR_2D_gau))=0;
    imagesc(FR_2D_gau,'AlphaData',imAlpha);
    Max_V                       = max(FR_2D_gau(:));
    Max_V                       = (floor(Max_V/2)+1)*2;
    caxis([0 Max_V]);  
    axis equal;
    xlim([1 40])
    ylim([1 40])
    axis tight
    axis off
    set(gcf,'color','w');
    %     c6                          = colorbar;
    %     title(c6,'R','FontSize',ClrBar_Font);
    %     title(['FR Htmp: Info=' num2str(info(segs))]);%
    hold on;
    FR_2D_gau                      = Spk_2D_pad_gau_depad_temp./occ_2D_pad_gau_depad;
    FR_2D_gau(FR_2D_gau==Inf)      = 0;
    FR_2D_gau_norm                 = FR_2D_gau/max(FR_2D_gau(:));
    FR_2D_gau_Binarized           = (FR_2D_gau_norm>=0.05);
    [B,L] = bwboundaries(FR_2D_gau_Binarized);
   for k=1:length(B)
   b = B{k};
   plot(b(:,2),b(:,1),'r','LineWidth',2);
   hold on;
   end 
   hold on;
end
drawnow;
imagewd = getframe(gcf);
imwrite(imagewd.cdata, [Folder 'Cells_Pos_Heatmap_Placefield.tif'],'WriteMode','append');
%%-------------------------------------------------------------------------
vw = VideoWriter(['D:\Cscope Data\PP367\Session1_12_51_06\Processed\' 'Occ_Neural_3_info_p_shf1st.avi'],'Motion JPEG AVI');
vw.Quality = 95;
open(vw)
tiff_info = imfinfo([Folder 'Occ_Neural_3.tif']);
 for ii             = 1 : size(tiff_info, 1)
    FR_2D_Norm_rgb  = double(imread([Folder 'Occ_Neural_3.tif'], ii));
    gcf=figure(700)
    imshow(uint8(FR_2D_Norm_rgb));
   
    title(['Info Content:' sprintf('%8.4f',info_Shuffle(1,ii)) ';P:' sprintf('%8.4f',infoP(ii))])
    disp([num2str(ii) '/' num2str(size(tiff_info, 1))]);
    imagewd = getframe(gcf);
    writeVideo(vw,imagewd);
end
close(vw);
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
load([Folder 'Info_P_Shf2nd.mat']);
vw = VideoWriter(['D:\Cscope Data\PP367\Session1_12_51_06\Processed\' 'Heatmap_T0-3_0-4_0-5.avi'],'Motion JPEG AVI');
vw.Quality                         = 95;
open(vw)

Place_Field_T= zeros(size(Spk_2D_pad_gau_depad,1),size(Spk_2D_pad_gau_depad,2));
for segs                           = 1:ns
    if (infoP(segs)>=0.95)
        gcf=figure(8002)
        gcf.Position                   = [180 180 1200 250];
        ClrBar_Font                    = 10;
        ha = tight_subplot(1,6,[.01 .01],[.01 .01],[.01 .01])
        sgtitle(['Cell: #' num2str(segs)]);
                
        Spk_2D_pad_gau_depad_temp      = Spk_2D_pad_gau_depad(:,:,segs);
        FR_2D_gau                      = Spk_2D_pad_gau_depad_temp./occ_2D_pad_gau_depad;
        FR_2D_gau_ref                  = FR_2D_gau;
        FR_2D_gau_ref(FR_2D_gau_ref==Inf)   = nan;
        FR_2D_gau(FR_2D_gau==Inf)      = 0;
        FR_2D_gau_norm                 = FR_2D_gau/max(FR_2D_gau(:));
        FR_2D_gau_Binarized_3          = (FR_2D_gau_norm>=0.2);
        FR_2D_gau_Binarized_4          = (FR_2D_gau_norm>=0.4);
        FR_2D_gau_Binarized_5          = (FR_2D_gau_norm>=0.5);
        
        Place_Field_T = Place_Field_T+ FR_2D_gau_norm;
                
        axes(ha(1))
        cm3                            = jet(255);
        %           cm3                            = colormap(ha(1),'jet');
        pos_1                          = get(ha(1),'position');
        pos_1(1)                       = pos_1(1);
        subplot('Position',pos_1);
        plot(posx,posy,'k','LineWidth' ,0.5);
        xlim([-binRad binRad])
        ylim([-binRad binRad])
        xticks(-binRad:5:binRad)
        yticks(-binRad:5:binRad)
        axis equal
        axis tight
        hold on;
        posx_spk_spd                   = posx(spks_binary(segs,:)&Speed_r_th_idx);
        posy_spk_spd                   = posy(spks_binary(segs,:)&Speed_r_th_idx);
        spks_spk_spd                   = spks(segs,spks_binary(segs,:)&Speed_r_th_idx)';
        spk_color                      = floor(spks_spk_spd/max(spks_spk_spd(:))*254)+1;
        scatter(posx_spk_spd,posy_spk_spd, 20,cm3(spk_color,:) ,'filled');
        %         colormap('jet');
        c3                             = colorbar;
        title(c3,'P','FontSize',ClrBar_Font);
        axis off;
        set(gcf,'color','w');          
        
        axes(ha(2))
        imAlpha=ones(size(FR_2D_gau_ref));
        imAlpha(isnan(FR_2D_gau_ref))=0;
        imagesc(FR_2D_gau_norm,'AlphaData',imAlpha);
        Max_V                       = max(FR_2D_gau_norm(:));
        Max_V                       = (floor(Max_V/2)+1)*2;
        %         caxis([0 Max_V-1]);
        %         caxis([0 5]);
        %         cm1                            = colormap(ha(1), 'parula')
        c1                             = colorbar;
        title(c1,'S','FontSize',8);
        title('Heatmap');
        axis equal;
        xlim([1 40])
        ylim([1 40])
        axis tight
        axis off        
        
        axes(ha(3))
        imagesc(FR_2D_gau_Binarized_3);
        Max_V                       = max(FR_2D_gau_Binarized_3(:));
        Max_V                       = (floor(Max_V/2)+1)*2;
        %         caxis([0 Max_V-1]);
        %         caxis([0 5]);
        %         cm1                            = colormap(ha(2), 'parula')
        c1                             = colorbar;
        title(c1,'S','FontSize',ClrBar_Font);
        title(['Heatmap:' num2str(0.2)]);
        axis equal;
        xlim([1 40])
        ylim([1 40])
        axis tight
        axis off
        
        axes(ha(4))
        imagesc(FR_2D_gau_Binarized_4);
        Max_V                       = max(FR_2D_gau_Binarized_4(:));
        Max_V                       = (floor(Max_V/2)+1)*2;
        %         caxis([0 Max_V-1]);
        %         caxis([0 5]);
        %         cm1                            = colormap(ha(3), 'parula')
        c1                             = colorbar;
        title(c1,'S','FontSize',ClrBar_Font);
        title(['Heatmap:' num2str(0.4)]);
        axis equal;
        xlim([1 40])
        ylim([1 40])
        axis tight
        axis off        
        
        axes(ha(5))
        imagesc( FR_2D_gau_Binarized_5);
        imagesc(FR_2D_gau_Binarized_4);
        Max_V                       = max(FR_2D_gau_Binarized_5(:));
        Max_V                       = (floor(Max_V/2)+1)*2;
        %         caxis([0 Max_V-1]);
        %         caxis([0 5]);
        %         cm1                            = colormap(ha(4), 'parula')
        c1                             = colorbar;
        title(c1,'S','FontSize',ClrBar_Font);
        title(['Heatmap:' num2str(0.5)]);
        axis equal;
        xlim([1 40])
        ylim([1 40])
        axis tight
        axis off        
        
        axes(ha(6))
        imAlpha=ones(size(Place_Field_T));
        imAlpha(isnan(Place_Field_T))=0;
        imagesc(Place_Field_T,'AlphaData',imAlpha);
        Max_V                       = max(Place_Field_T(:));
        Max_V                       = (floor(Max_V/2)+1)*2;
        %         caxis([0 Max_V-1]);
        %         caxis([0 5]);
        %         cm1                            = colormap(ha(4), 'parula')
        c1                             = colorbar;
        title(c1,'S','FontSize',ClrBar_Font);
        title(['Sum of Heatmaps']);
        axis equal;
        xlim([1 40])
        ylim([1 40])
        axis tight
        axis off
        drawnow;
        imagewd = getframe(gcf);
        writeVideo(vw,imagewd);
        close(gcf);
        
    end
end
close(vw);
%--------------------------------------------------------------------------
Place_Field_T                          = zeros(size(Spk_2D_pad_gau_depad,1),size(Spk_2D_pad_gau_depad,2));
sum(infoP>=0.95)/575
for segs                               = 1:ns
    
    if (infoP(segs)>=0.95)
        Spk_2D_pad_gau_depad_temp      = Spk_2D_pad_gau_depad(:,:,segs);
        FR_2D_gau                      = Spk_2D_pad_gau_depad_temp./occ_2D_pad_gau_depad;
        FR_2D_gau(FR_2D_gau==Inf)      = 0;
        FR_2D_gau_norm                 = FR_2D_gau/max(FR_2D_gau(:));
        FR_2D_gau_Binarized(:,:,segs)  = (FR_2D_gau_norm>=0.5); 
        Place_Field_T                  = Place_Field_T + FR_2D_gau_Binarized(:,:,segs);
        segs        
    end
    
end
    
for segs=1:ns
    gcf=figure(9001);
    FR_2D_gau_binar_temp              = double(FR_2D_gau_Binarized(:,:,segs))*2;
%     FR_2D_gau_binar_temp(FR_2D_gau_binar_temp==0)   = nan;
    imAlpha                           = ones(size(FR_2D_gau_binar_temp));
    imAlpha(isnan(FR_2D_gau_binar_temp))     = 0;
    imagesc(FR_2D_gau_binar_temp,'AlphaData',imAlpha);
    colorbar;
    Max_V                           = max(FR_2D_gau_binar_temp(:));
    Max_V                           = (floor(Max_V/2)+1)*2;
    caxis([0 Max_V]);  
    axis equal;
    xlim([1 40])
    ylim([1 40])
    axis tight
    axis off
    set(gcf,'color','w');
    imagewd = getframe(gcf);
    imwrite(imagewd.cdata, [Folder 'Place_Fields.tif'],'WriteMode','append')
    segs
end

    Place_Field_T(Place_Field_T==0)   = nan;
    imAlpha                           = ones(size(Place_Field_T));
    imAlpha(isnan(Place_Field_T))     = 0;
    gcf                               = figure(9002)
    imagesc(Place_Field_T,'AlphaData',imAlpha);
    colorbar;
    Max_V                             = max(Place_Field_T(:));
    Max_V                             = (floor(Max_V/2)+1)*2;
    caxis([0 Max_V]);  
    axis equal;
    xlim([1 40])
    ylim([1 40])
    axis tight
    axis off
    set(gcf,'color','w');
    imagewd                          = getframe(gcf);
    imwrite(imagewd.cdata, [Folder 'Cells_Heatmaps_Fig7h.tif']);    
    
    Place_Field_T_P                  = Place_Field_T;
    Place_Field_T_P(isnan(Place_Field_T_P))   = 0;
    Place_Field_T_P                  = Place_Field_T_P/sum(infoP>=0.95)*100;
    Place_Field_T_P(Place_Field_T_P==0)=nan;
    imAlpha                          = ones(size(Place_Field_T_P));
    imAlpha(isnan(Place_Field_T_P))  = 0;
    gcf=figure(9003)
    imagesc(Place_Field_T_P,'AlphaData',imAlpha);
    colorbar;
    Max_V                           = max(Place_Field_T_P(:));
    Max_V                           = (floor(Max_V/2)+1)*2;
    caxis([0 Max_V]);  
    axis equal;
    xlim([1 40])
    ylim([1 40])
    axis tight
    axis off
    set(gcf,'color','w');
    imagewd = getframe(gcf);
    imwrite(imagewd.cdata, [Folder 'Cells_Heatmaps_Fig7h.tif']);
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
load([Folder 'A_neuron_1_3d_arry_g.mat']);

A_Cells_T = A_Neuron_g;
A_Cells_T                       = A_Cells_T./max(A_Cells_T,[],[1 2]); 
Threshold = 0.5;
A_Cells_T(A_Cells_T<=Threshold)=0; %% Add a threshold for the footprin
A_Cells_T(isnan(A_Cells_T))=0;
for segno=1:size(A_Cells_T,3)
figure(777)
imshow((A_Cells_T(:,:,segno)*1))
pause 
end  
A_Cells_T_sum = sum(A_Cells_T,3);
A_Cells_T_sum = A_Cells_T_sum./max(A_Cells_T_sum(:));
find(A_Cells_T_sum==nan)
figure(778)
imshow(A_Cells_T_sum)
A_x    = size(A_Cells_T,1);
A_y    = size(A_Cells_T,2);
A_z    = size(A_Cells_T,3);


temp = hsv(256);
K    = size(A_Cells_T,3);
col  = temp(randi(256, K,1), :); % randi(64,K,1) generats KX1 matrix pseudorandom integers between 1 and 64.



A_2D   = reshape( A_Cells_T , A_x*A_y , A_z );
A_2D   = single(A_2D);

A_2D_rgb        = zeros( size(A_2D,1),size(A_2D,2), 3 );
A_2D_rgb(:,:,1) = A_2D*diag(col(:,1));
A_2D_rgb(:,:,2) = A_2D*diag(col(:,2));
A_2D_rgb(:,:,3) = A_2D*diag(col(:,3));
A_2D_rgb_T    = squeeze(sum(A_2D_rgb,2));
A_2D_rgb_T_3D = reshape(A_2D_rgb_T,A_x,A_y,3);
A_2D_rgb_T_3D = rot90(A_2D_rgb_T_3D);
% A_2D_rgb_T_3D = fliplr(A_2D_rgb_T_3D);
A_2D_rgb_T_3D = flipud(A_2D_rgb_T_3D);
figure(401)
imshow((A_2D_rgb_T_3D))
imwrite(uint8(A_2D_rgb_T_3D*255),[Folder 'A_Color_Total_Rand_Trd' num2str(Threshold) '.tif']);

BK = imread([Folder 'BK_Frame1.tif']);
figure(401)
imshow(A_2D_rgb_T_3D+repmat(double(BK),1,1,3)/255*1)
BK_2 = (~((A_Cells_T_sum')>0)).* double(BK);
figure(402)
imshow(BK_2/255)
figure(403)
imshow(A_2D_rgb_T_3D*1.2+repmat(double(BK_2),1,1,3)/255*1.5)

imwrite(uint8((A_2D_rgb_T_3D*1.2+repmat(double(BK_2),1,1,3)/255*1.5)*255),[Folder 'A_Color_Total_Rand_Trd' num2str(Threshold) '_BK.tif']);