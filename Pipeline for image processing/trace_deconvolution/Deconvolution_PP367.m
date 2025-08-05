clear all;
close all;
clc;

Folder ='D:\Cscope Data\PP367\Session1_12_51_06\Processed\';
addpath([Folder 'Matlab Codes\CNMF_E-master\OASIS_matlab\'],'-begin');
addpath([Folder 'Matlab Codes\CNMF_E-master\OASIS_matlab\functions\'],'-begin');
addpath([Folder 'Matlab Codes\CNMF_E-master\OASIS_matlab\packages\'],'-begin');
addpath([Folder 'Matlab Codes\CNMF_E-master\OASIS_matlab\packages\oasis\'],'-begin');

Folder1      = Folder;
load([Folder1 'C_Raw.mat']);

d1           = designfilt('lowpassfir','PassbandFrequency',0.1,'StopbandFrequency',0.4);
C_Raw_g      = C_Raw;
parfor cellno         = 1:size(C_Raw_g,1)   
C_Raw_g_temp          = C_Raw_g(cellno,:);
C_Raw_g(cellno,:)     = C_Raw_g_temp;
C_Raw_g_flt(cellno,:) = filtfilt(d1,C_Raw_g(cellno,:));
end

save([Folder1 'C_Raw_Filtered_ClG.mat'] , 'C_Raw_g_flt', '-v7.3')


dt_m                  = 0.0454;
yshift                = 1;

figure(100);
for cellno            = 1:20
    
    C_raw_fn  = C_Raw_g_flt(cellno,:);
    C_rawn    = C_Raw_g(cellno,:);
    
    Max       = max([max(C_raw_fn(:)),max(C_rawn(:))]);
    C_raw_fn  = C_raw_fn/Max;
    C_rawn    = C_rawn/Max;
    
    plot(C_rawn(1:1:size(C_rawn,2))+(cellno-1)*yshift,'Color',[100/255 100/255 100/255]);
    hold on;
    plot(C_raw_fn(1:1:size(C_rawn,2))+(cellno-1)*yshift,'r','LineWidth',0.2);
    axis tight
    
    xlim([1 size(C_rawn,2)])
    xticks(1:1*1000:size(C_rawn,2))
    xticklabels( string(dt_m*1:dt_m*1000:dt_m*size(C_rawn,2)) )
    hold on;    
end

tau_dr                  = [45,0.8];
window                  = 200;
ht                      = exp2kernel(tau_dr, window);
figure(1)
plot(ht)
lambda                  = 8;
parfor segs             = 1:size(C_Raw_g_flt,1)
[c_Dnl, s_Dnl, options] = deconvolveCa(C_Raw_g_flt(segs,:), 'exp2', tau_dr, 'foopsi','optimize_pars',false,'optimize_b',true,'lambda',lambda,'maxIter',5,'window',300,'shift',100)
c_2d(segs,:)            = c_Dnl;
s_2d(segs,:)            = s_Dnl;
options_cells(segs,:)   = options;
segs
end

save([Folder1 'Spk_Daniel_g.mat'],'s_2d','options_cells','c_2d','-v7.3');

f=figure(2000);
sgtitle(['RisingT:' num2str(tau_dr(2)*1) 'frm; DecayT:' num2str(tau_dr(1)*1) 'frm']);
f.Position = [10 10 1000 1000];
yshift     = 1;

for cellno = 1:50
    
    C_raw_fn      = C_Raw_g_flt(cellno,:);
    C_rawn        = C_Raw_g(cellno,:);
    
    c_2dn         = c_2d(cellno,:);
    s_2dn         = s_2d(cellno,:);
    
    Max           = max([max(C_raw_fn(:)),max(C_rawn(:)),max(c_2dn(:)),max(s_2dn(:))]);
    C_raw_fn      = C_raw_fn/Max;
    C_rawn        = C_rawn/Max;
    
    c_2dn         = c_2dn/Max;
    s_2dn         = s_2dn/Max;
    
    plot(s_2dn(1:size(c_2dn,2))+(cellno-1)*yshift,'r','LineWidth',1);
    hold on;
    plot(c_2dn(1:size(c_2dn,2))+(cellno-1)*yshift,'k','LineWidth',0.2);
    hold on;
    
%     plot(C_rawn(1:1:size(c_2dn,2))+(cellno-1)*yshift,'Color',[100/255 100/255 100/255]);
    plot(C_raw_fn(1:1:size(c_2dn,2))+(cellno-1)*yshift,'b','LineWidth',0.5);
    axis tight
    
    xlim([1 size(c_2dn,2)])
    xticks(1:1*1000:size(c_2dn,2))
    xticklabels( string(dt_m*1:dt_m*1000:dt_m*size(c_2dn,2)) )
    hold on;
    
end
