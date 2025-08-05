% clear all;
% close all;
% clc;

% Folder ='F:\2-Cscope Project-FOV-3mm_NA-0d1\Cscope Exp\Mice_New_LinearTrack\PP367_Ses1_11_09_30\Processed\';
addpath(genpath( 'OASIS_matlab'), '-begin');

% load([Folder 'A_T_clean.mat']);

d1           = designfilt('lowpassfir','PassbandFrequency',0.3,'StopbandFrequency',0.4 );
C_Raw_g      = detrended_trace(A_neuron_good_idx+1,:);
parfor cellno         = 1:size(C_Raw_g,1)   
C_Raw_g_flt(cellno,:) = filtfilt(d1,C_Raw_g(cellno,:));
end


dt_m                  = 1;
yshift                = 1;

idrange=5050:5060;

figure(100);
% for cellno            = idrange
% 
%     C_raw_fn  = C_Raw_g_flt(cellno,:);
%     C_rawn    = C_Raw_g(cellno,:);
% 
%     Max       = max([max(C_raw_fn(:)),max(C_rawn(:))]);
%     C_raw_fn  = C_raw_fn/Max;
%     C_rawn    = C_rawn/Max;
% 
%     plot(C_rawn(1:1:size(C_rawn,2))+(cellno-1)*yshift,'Color',[100/255 100/255 100/255]);
%     hold on;
%     plot(C_raw_fn(1:1:size(C_rawn,2))+(cellno-1)*yshift,'r','LineWidth',0.2);
%     axis tight
% 
%     xlim([1 size(C_rawn,2)])
%     xticks(1:1*1000:size(C_rawn,2))
%     xticklabels( string(dt_m*1:dt_m*1000:dt_m*size(C_rawn,2)) )
%     hold on;    
% end

tau_dr                  = [10,0];
window                  = 200;
ht                      = exp2kernel(tau_dr, window);
figure(1)
plot(ht)
lambda                  =6;

parfor segs             = 1:size(C_Raw_g_flt,1)
% [c_Dnl, s_Dnl, options] = deconvolveCa(C_Raw_g_flt(segs,:), 'exp2', tau_dr, 'foopsi','optimize_pars',true,'optimize_b',true,'lambda',lambda,'maxIter',10,'window',300,'shift',100)
[c_Dnl, s_Dnl, options] = deconvolveCa(C_Raw_g_flt(segs,:), 'exp2',tau_dr,'foopsi','optimize_pars',0,'optimize_b',true,'lambda',lambda,'maxIter',10,'window',300,'shift',100)
c_2d(segs,:)            = c_Dnl;
s_2d(segs,:)            = s_Dnl;
options_cells(segs,:)   = options;
segs
end

% save([Folder 'A_C_S_T_clean.mat'] ,'A_good_3D','A_bad_list','C_Raw_clean', 's_2d','options_cells','c_2d','-v7.3')
save([Folder 'A_C_S_T_clean.mat'] ,'A_neuron_bad_idx', 'A_neuron_good_idx', 'A_neuron_sparse', 'coordinates', 'C_trace', 'C_raw','detrended_trace', 's_2d','options_cells','c_2d','-v6')
f=figure(2000);
sgtitle(['RisingT:' num2str(tau_dr(2)*1) 'frm; DecayT:' num2str(tau_dr(1)*1) 'frm']);
f.Position = [10 10 1000 1000];

for cellno = idrange
    
    C_raw_fn      = C_Raw_g_flt(cellno,:);
    C_rawn        = C_Raw_g(cellno,:);
    
    c_2dn         = c_2d(cellno,:);
    s_2dn         = s_2d(cellno,:);
    
    % Max           = max([max(C_raw_fn(:)),max(C_rawn(:)),max(c_2dn(:)),max(s_2dn(:))]);
    % C_raw_fn      = C_raw_fn/Max;
    % C_rawn        = C_rawn/Max;
    % 
    % c_2dn         = c_2dn/Max;
    % s_2dn         = s_2dn/Max;
    Max           = max([max(C_raw_fn(:)),max(C_rawn(:)),max(c_2dn(:))]);
    C_raw_fn      = C_raw_fn/Max;
    C_rawn        = C_rawn/Max;

    c_2dn         = c_2dn/Max;
    s_2dn         = s_2dn/max(s_2dn(:));
    
    plot(s_2dn(1:size(c_2dn,2))+(cellno-1)*yshift,'r','LineWidth',1);
    hold on;
    plot(c_2dn(1:size(c_2dn,2))+(cellno-1)*yshift,'k','LineWidth',0.2);
    hold on;

%     plot(C_rawn(1:1:size(c_2dn,2))+(cellno-1)*yshift,'Color',[100/255 100/255 100/255]);
    plot(C_raw_fn(1:1:size(c_2dn,2))+(cellno-1)*yshift,'b','LineWidth',0.5);
    axis tight
    
    xlim([1 size(c_2dn,2)])
    xticks(1:1*100:size(c_2dn,2))
    xticklabels( string(dt_m*1:dt_m*100:dt_m*size(c_2dn,2)) )
    hold on;
    
end

% f=figure(2001);
% sgtitle(['RisingT:' num2str(tau_dr(2)*1) 'frm; DecayT:' num2str(tau_dr(1)*1) 'frm']);
% f.Position = [10 10 1000 1000];
% yshift     = 1;
% 
% for cellno = idrange
% 
%     C_raw_fn      = C_Raw_g_flt(cellno,:);
%     C_rawn        = C_Raw_g(cellno,:);
% 
%     c_2dn         = c_2d(cellno,:);
%     s_2dn         = s_2d(cellno,:);
% 
%     Max           = max([max(C_raw_fn(:)),max(C_rawn(:)),max(c_2dn(:)),max(s_2dn(:))]);
%     C_raw_fn      = C_raw_fn/Max;
%     C_rawn        = C_rawn/Max;
% 
%     c_2dn         = c_2dn/Max;
%     s_2dn         = s_2dn/Max;
% 
%     plot(C_raw_fn(1:size(c_2dn,2))-c_2dn(1:1:size(c_2dn,2))+(cellno-1)*yshift,'r','LineWidth',1);
%     hold on;
%     axis tight
% 
%     xlim([1 size(c_2dn,2)])
%     xticks(1:1*1000:size(c_2dn,2))
%     xticklabels( string(dt_m*1:dt_m*1000:dt_m*size(c_2dn,2)) )
%     hold on;
% 
% end

% figure
% plot([options_cells(idrange).sn])
