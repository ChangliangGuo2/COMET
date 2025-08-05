%%=========================================================================
%%=========================================================================
%%
%%  ---- This code is for analyzing the fear conditionning experiment 
%%
%%  ---- by Changliang Guo 
%%
%%  ---- 05/15/2024
%%
%%=========================================================================
%%=========================================================================

clc;
clear;
close all;

%%-------------------------------------------------------------------------

%% Load the Mouse_STRUCT of 5 mice to one STRUCT

for mouseidx = 1:5

    Dir                           = 'D:\XM\';
    
    MiceName_T                    = {'FC_M94' ;'FC_M367'; 'FC_M369'; 'FC_M370' ;'FC_M371'};
    MiceName                      = MiceName_T{mouseidx};
    disp(['Processing Mouse: ' MiceName])
    
    path                          = [Dir MiceName '\'];
    
    filenames                     = [path MiceName '_Mouse_STRUCT.mat'];


    load(  filenames  );
    Mice_STRUCT_T(mouseidx) = Mouse_STRUCT;

    mouseidx

end

imwrite(uint8(Mice_STRUCT_T(2).ATLAS90*255) , "D:\XM\FC_M367\ATLAS.tiff")

clear  Mouse_STRUCT ;

save('D:\XM\FC_M94\Matlab codes\Mice_STRUCT_T.mat',"Mice_STRUCT_T","-v7.3");
%%-------------------------------------------------------------------------

Regions_overlap = intersect( Mice_STRUCT_T(1).REGION_NAME(Mice_STRUCT_T( 1 ).REGION_ID>0),  Mice_STRUCT_T(1).REGION_NAME(Mice_STRUCT_T( 1 ).REGION_ID>0));

for mouseidx=1:4

    Regions_overlap = intersect(Regions_overlap,  Mice_STRUCT_T(mouseidx).REGION_NAME(Mice_STRUCT_T( mouseidx ).REGION_ID>0));

end


Mice_STRUCT_T_Regions_overlap_Num_LR = zeros(14,1);
Mice_STRUCT_T_Regions_overlap_per_LR = zeros(14,1);

for mouseidx=1:5

    Regions_overlap_index{mouseidx} = ismember(Mice_STRUCT_T(mouseidx).REGION_NAME , Regions_overlap);
    
    Regions_overlap_inseq           = Mice_STRUCT_T(mouseidx).REGION_NAME(Regions_overlap_index{mouseidx} )
    
    Mice_STRUCT_T_Regions_overlap_Num(:,mouseidx) = Mice_STRUCT_T( mouseidx ).REGION_ID(Regions_overlap_index{mouseidx});
    Mice_STRUCT_T_Regions_overlap_per(:,mouseidx) = Mice_STRUCT_T_Regions_overlap_Num(:,mouseidx)./size(Mice_STRUCT_T(mouseidx).Neuron_REGION_ID,1);
    
    Mice_STRUCT_T_Regions_overlap_Num_LR(1:13,mouseidx) = Mice_STRUCT_T_Regions_overlap_Num(1:2:25,mouseidx) + Mice_STRUCT_T_Regions_overlap_Num(2:2:26,mouseidx);
    Mice_STRUCT_T_Regions_overlap_per_LR(1:13,mouseidx) = Mice_STRUCT_T_Regions_overlap_per(1:2:25,mouseidx) + Mice_STRUCT_T_Regions_overlap_per(2:2:26,mouseidx);
    
    Mice_STRUCT_T_Regions_overlap_Num_LR(14,mouseidx) = Mice_STRUCT_T_Regions_overlap_Num(27,mouseidx) ;
    Mice_STRUCT_T_Regions_overlap_per_LR(14,mouseidx) = Mice_STRUCT_T_Regions_overlap_per(27,mouseidx) ;

end

Regions_overlap_LR = Regions_overlap_inseq(1:2:end)


Mice_STRUCT_T_Regions_overlap_per_LR_mean = mean(Mice_STRUCT_T_Regions_overlap_per_LR,2)

[value, index] = sort(Mice_STRUCT_T_Regions_overlap_per_LR_mean,'descend')

Mice_STRUCT_T_Regions_overlap_per_LR_descend = Mice_STRUCT_T_Regions_overlap_per_LR(index,:);
Mice_STRUCT_T_Regions_overlap_Num_LR_descend = Mice_STRUCT_T_Regions_overlap_Num_LR(index,:);


Regions_overlap_LR_descend = Regions_overlap_LR(index)


color_inuse = [
    
   
    255,99,71;

   
    124,252,0;
    
  
    0,206,209;

    127,255,212;
   

    220,20,60;
  
    
    102,205,170;
   

    250,128,114;
   

    233,150,122;
  

    138,43,226;
  

    100,149,237;
   

    32,178,170;
   
 
    188,143,143;
  

    123,104,238;
    
    255,105,180;
  
  
    ]/255;


%%-------------------------------------------------------------------------
%% Chose a mouse for analysis first


MOUSE_CHOSEN                   = 1

Neuron_REGION_ID_index         = Mice_STRUCT_T( MOUSE_CHOSEN ).Neuron_REGION_ID;
[Bvalue,Neuron_REGION_ID_Sort] = sort(Neuron_REGION_ID_index);

Bvalue_unique                  = unique(Bvalue)
[Bvalue_unique_NUM n]          = hist(Bvalue,Bvalue_unique)

load(['D:\XM\FC_M94\Matlab codes\slanCM\' 'slanCM_Data.mat'])

Mice_STRUCT_T( MOUSE_CHOSEN ).REGION_NAME(  Mice_STRUCT_T(MOUSE_CHOSEN).REGION_ID>0  )

Trace_MOUSE                    = Mice_STRUCT_T( MOUSE_CHOSEN ).NEURON.Trace_unify(Neuron_REGION_ID_Sort,:);

close all


%-----------------------------------------


f1                    = figure(1)
hFig                  = gcf;
hAx                   = gca;

Trace_MOUSE_2         = imadjust(Trace_MOUSE, [0.05 1], [0 0.7], 1);
imtrace               = imshow((Trace_MOUSE_2*5))
colormap( slandarerCM(1).Colors{5} )

xLimits               = xlim;                                % 获取 x 轴限制
yLimits               = ylim;

Trace_unify_mean      = mean( Mice_STRUCT_T( MOUSE_CHOSEN ).NEURON.Trace_unify);
Trace_unify_mean_norm = rescale(Trace_unify_mean);


%-----------------------------------------


f2                    = figure('Position', get(f1, 'Position')); % 确保大小和位置一致
set(gca, 'Position', get(hAx, 'Position'))
% %
CueShock_marker_color(Mice_STRUCT_T( MOUSE_CHOSEN ).NEURON.CUE_Frame(1,:), Mice_STRUCT_T( MOUSE_CHOSEN ).NEURON.CUE_Frame(2,:),Mice_STRUCT_T( MOUSE_CHOSEN ).NEURON.SHOCK_Frame(1,:),Mice_STRUCT_T( MOUSE_CHOSEN ).NEURON.SHOCK_Frame(2,:),size(Mice_STRUCT_T( MOUSE_CHOSEN ).Neuron_REGION_ID,1));
% %
hold on
plot(  Mice_STRUCT_T( MOUSE_CHOSEN ).NEURON.Speed*50,'LineWidth',1,'Color','k'   )
hold on
plot(  Trace_unify_mean_norm*1000+1100,'LineWidth',1,'Color','red'  )
hold on
yline(1100,'r--','LineWidth',1.5 )

xlim(xLimits);                                
ylim(yLimits);                                 
set(gca, 'XTickLabel', []); 
set(gca, 'YTickLabel', []);  
xlabel('');  
ylabel('');  
set(gca, 'XColor', 'none') 
set(gca, 'YColor', 'none')
set(gca, 'Box', 'off');  


%-----------------------------------------


f3 = figure('Position', get(f1, 'Position')); % 确保大小和位置一致
set(gca, 'Position', get(hAx, 'Position'))
% %
CueShock_marker_color(Mice_STRUCT_T(MOUSE_CHOSEN).NEURON.CUE_Frame(1,:), Mice_STRUCT_T(MOUSE_CHOSEN).NEURON.CUE_Frame(2,:),Mice_STRUCT_T(MOUSE_CHOSEN).NEURON.SHOCK_Frame(1,:),Mice_STRUCT_T(MOUSE_CHOSEN).NEURON.SHOCK_Frame(2,:),size(Mice_STRUCT_T(MOUSE_CHOSEN).Neuron_REGION_ID,1));
% %
hold on
plot(  Trace_unify_mean_norm*1000,'LineWidth',1,'Color','red'  )
xlim(xLimits);                               
ylim(yLimits);                              
set(gca, 'XTickLabel', []); 
set(gca, 'YTickLabel', []);  
xlabel(''); 
ylabel('');  
set(gca, 'XColor', 'none') 
set(gca, 'YColor', 'none')
set(gca, 'Box', 'off'); 


%-----------------------------------------


Bvalue_unique_NUM_CUMSUM = cumsum(Bvalue_unique_NUM);
Bvalue_unique_NUM        = Bvalue_unique_NUM'
Bvalue_unique_NAME       = {Mice_STRUCT_T( MOUSE_CHOSEN ).REGION_NAME{Bvalue_unique}  }


%-----------------------------------------


f4 = figure('Position', get(f1, 'Position')); % 确保大小和位置一致
set(gca, 'Position', get(hAx, 'Position'))

%     hFig = gcf;
%     hAx  = gca;
%     % set the figure to full screen
%     set(hFig,'units','normalized','outerposition',[0 0 0.5 1]);
%     % set the axes to full screen
%     set(hAx,'Unit','normalized','Position',[0 0 0.5 1]);
% xline(1)

xticks(0:1)
ylim([1 Bvalue_unique_NUM_CUMSUM(end)])
yticks([1,Bvalue_unique_NUM_CUMSUM])

% yticklabels(string(Bvalue_unique_NUM_CUMSUM))


%--------------------------------------------------------------------------



figure (5)
CueShock_marker_color(Mice_STRUCT_T(MOUSE_CHOSEN).NEURON.CUE_Frame(1,:), Mice_STRUCT_T(MOUSE_CHOSEN).NEURON.CUE_Frame(2,:),Mice_STRUCT_T(MOUSE_CHOSEN).NEURON.SHOCK_Frame(1,:),Mice_STRUCT_T(MOUSE_CHOSEN).NEURON.SHOCK_Frame(2,:),1);
% %
hold on
plot(  Trace_unify_mean_norm*1,'LineWidth',1,'Color','red'  )
xlim(xLimits);                                                             
set(gca, 'XTickLabel', []); 
set(gca, 'YTickLabel', []);  
xlabel(''); 
ylabel('');  
set(gca, 'XColor', 'none') 
set(gca, 'YColor', 'none')
set(gca, 'Box', 'off'); 



%--------------------------------------------------------------------------

%%=========================================================================


%%----------List the cortex regions by M,SS,DSP,VIS------------------------

Regions_GROUP = {
     'MOs1_R'   ,
     'MOs1_L'   ,
     'MOp1_R'   ,
     'MOp1_L'   ,
    'SSp-ul1_R' ,
    'SSp-ul1_L' ,
    'SSp-ll1_R' ,
    'SSp-ll1_L' ,
    'SSp-un1_R' ,
    'SSp-un1_L' ,
    'SSp-tr1_R' ,
    'SSp-tr1_L' ,
    'SSp-bfd1_R',
    'SSp-bfd1_L',
    'RSPd1_R'   ,
    'RSPd1_L'   ,
    'RSPagl1_R' ,
    'RSPagl1_L' ,
    'RSPv1_L'   ,
    'VISa1_R'   ,
    'VISa1_L'   ,  
    'VISrl1_R'  ,
    'VISrl1_L'  ,
    'VISam1_R'  ,
    'VISam1_L'  ,
    'VISpm1_R'  ,
    'VISpm1_L'  ,  
    'VISp1_R'   ,
    'VISp1_L'   ,
    'VISal1_L'  };
  

for i = 1:length(Regions_GROUP)

    Regions_GROUP_Index(i) = find(strcmp(Regions_GROUP(i), Mice_STRUCT_T( MOUSE_CHOSEN ).REGION_NAME));

end


Trace_M94_ReGroup      = [];
RawAct_M94_ReGroup     = [];


for i                  = 1:length(Regions_GROUP_Index)

    Trace_M94_ReGroup  = [Trace_M94_ReGroup; Mice_STRUCT_T(MOUSE_CHOSEN).NEURON.Trace_unify( Mice_STRUCT_T(MOUSE_CHOSEN).Neuron_REGION_ID==Regions_GROUP_Index(i) , : )];
    
    RawAct_M94_ReGroup = [RawAct_M94_ReGroup;Mice_STRUCT_T(MOUSE_CHOSEN).NEURON.RawAct_unify( Mice_STRUCT_T(MOUSE_CHOSEN).Neuron_REGION_ID==Regions_GROUP_Index(i) , : )];

end


figure(10)
    plot(mean(RawAct_M94_ReGroup))
    hold on
    plot(mean(Mice_STRUCT_T(MOUSE_CHOSEN).NEURON.RawAct_unify))


[~, index] = ismember(Regions_GROUP_Index, Bvalue_unique);


%%=========================================================================

close all
f1   = figure(20)
hFig = gcf;
hAx  = gca;

%     set(hAx,'position',[100,100,500,1000])    
%     set the figure to full screen
%     set(hFig,'units','normalized','outerposition',[0 0 0.6 1]);
%     % set the axes to full screen
%     set(hAx,'Unit','normalized','Position',[0 0 0.6 1]);

Trace_M94_ReGroup_2 = imadjust(Trace_M94_ReGroup, [0.05 1], [0 0.7], 1);
imshow( (Trace_M94_ReGroup_2*5)  )
colormap(  slandarerCM(1).Colors{5} )

%     hold on
%     xline(Mice_STRUCT_T(5).NEURON.CUE_Frame(1,:),'Color',[1 1 1],'LineWidth',1     )
%     hold on
%     xline(Mice_STRUCT_T(5).NEURON.CUE_Frame(2,:),'Color',[1 1 1],'LineWidth',1     )
%     hold on
%     xline(Mice_STRUCT_T(5).NEURON.SHOCK_Frame(1,:),'Color',[1 1 1],'LineWidth',1   )
%     hold on
%     xline(Mice_STRUCT_T(5).NEURON.SHOCK_Frame(2,:),'Color',[1 1 1],'LineWidth',0.5 )

xLimits               = xlim;                                % 获取 x 轴限制
yLimits               = ylim;

Trace_unify_mean      = mean( Mice_STRUCT_T(MOUSE_CHOSEN).NEURON.Trace_unify);
Trace_unify_mean_norm = rescale(Trace_unify_mean);



f2                    = figure('Position', get(f1, 'Position')); % 确保大小和位置一致

set(gca, 'Position', get(hAx, 'Position'))

% %--------------------------------------------------------------------------
CueShock_marker_color(Mice_STRUCT_T(MOUSE_CHOSEN).NEURON.CUE_Frame(1,:), Mice_STRUCT_T(MOUSE_CHOSEN).NEURON.CUE_Frame(2,:),Mice_STRUCT_T(MOUSE_CHOSEN).NEURON.SHOCK_Frame(1,:),Mice_STRUCT_T(MOUSE_CHOSEN).NEURON.SHOCK_Frame(2,:),size(Mice_STRUCT_T(MOUSE_CHOSEN).Neuron_REGION_ID,1));
% %--------------------------------------------------------------------------
hold on
plot(  Mice_STRUCT_T(MOUSE_CHOSEN).NEURON.Speed*50,'LineWidth',1,'Color','k'  )
hold on
plot(  Trace_unify_mean_norm*1000+1100,'LineWidth',1,'Color','red'  )
hold on
yline(1100,'r--','LineWidth',1.5 )

xlim(xLimits);                                 % 设置 x 轴限制
ylim(yLimits);                                 % 设置 y 轴限制

set(gca, 'XTickLabel', []);  % 去掉 x 轴的刻度标签

set(gca, 'YTickLabel', []);  % 去掉 y 轴的刻度标签

% 去掉 x 和 y 轴的标签

xlabel('');  % 去掉 x 轴标签

ylabel('');  % 去掉 y 轴标签
set(gca, 'XColor', 'none') 
set(gca, 'YColor', 'none')
set(gca, 'Box', 'off');  % 关闭边框



f3 = figure('Position', get(f1, 'Position')  ); % 确保大小和位置一致
set(gca,    'Position', get(hAx, 'Position') )
% %--------------------------------------------------------------------------
CueShock_marker_color(Mice_STRUCT_T(MOUSE_CHOSEN).NEURON.CUE_Frame(1,:), Mice_STRUCT_T(MOUSE_CHOSEN).NEURON.CUE_Frame(2,:),Mice_STRUCT_T(MOUSE_CHOSEN).NEURON.SHOCK_Frame(1,:),Mice_STRUCT_T(MOUSE_CHOSEN).NEURON.SHOCK_Frame(2,:),size(Mice_STRUCT_T(MOUSE_CHOSEN).Neuron_REGION_ID,1));
% %--------------------------------------------------------------------------
hold on
plot(  Trace_unify_mean_norm*1000,'LineWidth',1,'Color','red'  )
xlim(xLimits);                                 % 设置 x 轴限制
ylim(yLimits);                                 % 设置 y 轴限制
set(gca, 'XTickLabel', []);  % 去掉 x 轴的刻度标签

set(gca, 'YTickLabel', []);  % 去掉 y 轴的刻度标签

% 去掉 x 和 y 轴的标签

xlabel('');  % 去掉 x 轴标签

ylabel('');  % 去掉 y 轴标签
set(gca, 'XColor', 'none') 
set(gca, 'YColor', 'none')
set(gca, 'Box', 'off');  % 关闭边框


Bvalue_unique_NUM_ReGroup        = Bvalue_unique_NUM(index(find(index)))

Bvalue_unique_NUM_ReGroup_CUMSUM = cumsum(Bvalue_unique_NUM_ReGroup)';
Bvalue_unique_NUM_ReGroup        = Bvalue_unique_NUM_ReGroup'



f4 = figure('Position', get(f1, 'Position') ); % 确保大小和位置一致
set(gca,    'Position', get(hAx, 'Position'))

%     hFig = gcf;
%     hAx  = gca;
%     % set the figure to full screen
%     set(hFig,'units','normalized','outerposition',[0 0 0.5 1]);
%     % set the axes to full screen
%     set(hAx,'Unit','normalized','Position',[0 0 0.5 1]);
%     xline(1)

xticks( 0:1 )
ylim( [1 Bvalue_unique_NUM_ReGroup_CUMSUM(end)] )
yticks( [1,Bvalue_unique_NUM_ReGroup_CUMSUM]    )
yticklabels(string(Bvalue_unique_NUM_ReGroup)   )

Bvalue_unique_NAME    = { Mice_STRUCT_T(MOUSE_CHOSEN).REGION_NAME{Bvalue_unique}  }

%==========================================================================
%==========================================================================

%Show MO., SSP., RSP. and VIS. 


f100 = figure(100); % 确保大小和位置一致
f100.Position = [100 100 800 200];
% %
CueShock_marker_color(Mice_STRUCT_T(MOUSE_CHOSEN).NEURON.CUE_Frame(1,:), Mice_STRUCT_T(MOUSE_CHOSEN).NEURON.CUE_Frame(2,:),Mice_STRUCT_T(MOUSE_CHOSEN).NEURON.SHOCK_Frame(1,:),Mice_STRUCT_T(MOUSE_CHOSEN).NEURON.SHOCK_Frame(2,:),size(Mice_STRUCT_T(MOUSE_CHOSEN).Neuron_REGION_ID,1));
% %
hold on

plot(  Trace_unify_mean*1,'LineWidth',1,'Color','k'  )
hold on
Trace_M94_ReGroup_mean = (mean(Trace_M94_ReGroup(1:Bvalue_unique_NUM_ReGroup_CUMSUM(4),:)));

plot( Trace_M94_ReGroup_mean*1,'LineWidth',1.5,'Color',[9 153 99]/255  )
hold on;

% RawAct_M94_ReGroup_mean = (mean(RawAct_M94_ReGroup(1:Bvalue_unique_NUM_ReGroup_CUMSUM(4),:)));
% plot( RawAct_M94_ReGroup_mean*1,'LineWidth',1,'Color',[255 0 0]/255*1  )
% 
% hold on
% plot(  Mice_STRUCT_T(MOUSE_CHOSEN).NEURON.RawAct_unify_mean*1,'LineWidth',1,'Color','k'  )

xlim( xLimits);                                 % 设置 x 轴限制
ylim( [0,0.2]  );                                 % 设置 y 轴限制
set(gca, 'XTickLabel', []);  % 去掉 x 轴的刻度标签
set(gca, 'YTickLabel', []);  % 去掉 y 轴的刻度标签
% 去掉 x 和 y 轴的标签
xlabel('');  % 去掉 x 轴标签
ylabel('');  % 去掉 y 轴标签
set(gca, 'XColor', 'none') 
set(gca, 'YColor', 'none')
set(gca, 'Box', 'off');  % 关闭边框

%--------------------------------------------------------------------------

f101 = figure(101); % 确保大小和位置一致
f101.Position = [100 100 800 200];
% %
CueShock_marker_color(Mice_STRUCT_T(5).NEURON.CUE_Frame(1,:), Mice_STRUCT_T(5).NEURON.CUE_Frame(2,:),Mice_STRUCT_T(5).NEURON.SHOCK_Frame(1,:),Mice_STRUCT_T(5).NEURON.SHOCK_Frame(2,:),size(Mice_STRUCT_T(5).Neuron_REGION_ID,1));
% %
hold on
plot(  Trace_unify_mean*1,'LineWidth',1,'Color','k'  )
hold on
Trace_M94_ReGroup_mean = (mean(Trace_M94_ReGroup(Bvalue_unique_NUM_ReGroup_CUMSUM(4)+1:Bvalue_unique_NUM_ReGroup_CUMSUM(14),:)));
plot( Trace_M94_ReGroup_mean*1,'LineWidth',1.5,'Color',[85 160 251]/255  )
xlim(xLimits);                                 % 设置 x 轴限制
ylim([0,0.2]);                                 % 设置 y 轴限制
set(gca, 'XTickLabel', []);  % 去掉 x 轴的刻度标签
set(gca, 'YTickLabel', []);  % 去掉 y 轴的刻度标签
% 去掉 x 和 y 轴的标签
xlabel('');  % 去掉 x 轴标签
ylabel('');  % 去掉 y 轴标签
set(gca, 'XColor', 'none') 
set(gca, 'YColor', 'none')
set(gca, 'Box', 'off');  % 关闭边框


%--------------------------------------------------------------------------

f102 = figure(102); % 确保大小和位置一致
f102.Position = [100 100 800 200];
% %
CueShock_marker_color(Mice_STRUCT_T(5).NEURON.CUE_Frame(1,:), Mice_STRUCT_T(5).NEURON.CUE_Frame(2,:),Mice_STRUCT_T(5).NEURON.SHOCK_Frame(1,:),Mice_STRUCT_T(5).NEURON.SHOCK_Frame(2,:),size(Mice_STRUCT_T(5).Neuron_REGION_ID,1));
% %
hold on
plot(  Trace_unify_mean*1,'LineWidth',1,'Color','k'  )
hold on
Trace_M94_ReGroup_mean = (mean(Trace_M94_ReGroup(Bvalue_unique_NUM_ReGroup_CUMSUM(14)+1:Bvalue_unique_NUM_ReGroup_CUMSUM(19),:)));
plot( Trace_M94_ReGroup_mean*1,'LineWidth',1.5,'Color',[255 128 0]/255  )
xlim(xLimits);                                 % 设置 x 轴限制
ylim([0,.2]);                                 % 设置 y 轴限制
set(gca, 'XTickLabel', []);  % 去掉 x 轴的刻度标签
set(gca, 'YTickLabel', []);  % 去掉 y 轴的刻度标签
% 去掉 x 和 y 轴的标签
xlabel('');  % 去掉 x 轴标签
ylabel('');  % 去掉 y 轴标签
set(gca, 'XColor', 'none') 
set(gca, 'YColor', 'none')
set(gca, 'Box', 'off');  % 关闭边框

%--------------------------------------------------------------------------

f103 = figure(103); % 确保大小和位置一致
f103.Position = [100 100 800 200];
% %
CueShock_marker_color(Mice_STRUCT_T(5).NEURON.CUE_Frame(1,:), Mice_STRUCT_T(5).NEURON.CUE_Frame(2,:),Mice_STRUCT_T(5).NEURON.SHOCK_Frame(1,:),Mice_STRUCT_T(5).NEURON.SHOCK_Frame(2,:),size(Mice_STRUCT_T(5).Neuron_REGION_ID,1));
% %
hold on
plot(  Trace_unify_mean*1,'LineWidth',1,'Color','k'  )
hold on
Trace_M94_ReGroup_mean = (mean(Trace_M94_ReGroup(Bvalue_unique_NUM_ReGroup_CUMSUM(25)+1:end,:)));
plot( Trace_M94_ReGroup_mean*1,'LineWidth',1.5,'Color',[249 100 149]/255  )

% RawAct_M94_ReGroup_mean = (mean(RawAct_M94_ReGroup(Bvalue_unique_NUM_ReGroup_CUMSUM(25)+1:end,:)));
% plot( RawAct_M94_ReGroup_mean*1,'LineWidth',1,'Color',[255 0 0]/255*1  )

% hold on
% plot(   Mice_STRUCT_T(5).NEURON.RawAct_unify_mean*1,'LineWidth',1,'Color','k'  )

xlim(xLimits);                                 % 设置 x 轴限制
ylim([0,.2]);                                 % 设置 y 轴限制
set(gca, 'XTickLabel', []);  % 去掉 x 轴的刻度标签
set(gca, 'YTickLabel', []);  % 去掉 y 轴的刻度标签
% 去掉 x 和 y 轴的标签
xlabel('');  % 去掉 x 轴标签
ylabel('');  % 去掉 y 轴标签
set(gca, 'XColor', 'none') 
set(gca, 'YColor', 'none')
set(gca, 'Box', 'off');  % 关闭边框

%==========================================================================
%==========================================================================

REGIOIN_EXSIT         = (~isnan(Mice_STRUCT_T(1).NEURON.Trace_unify_Region_Avg(:,1)));
[Cortex_Incommon]     = Mice_STRUCT_T(1).REGION_NAME(REGIOIN_EXSIT);

for mouseidx          = 1:4

    REGIOIN_EXSIT     = (~isnan(Mice_STRUCT_T(mouseidx+1).NEURON.Trace_unify_Region_Avg(:,1)));
    [Cortex_Incommon] = intersect(Cortex_Incommon, Mice_STRUCT_T(mouseidx+1).REGION_NAME(REGIOIN_EXSIT)) ;

end

CortexID_Mouse                    = []

for mouseidx                      = 1:5
    
    [Cortex_Incommon,ai,CortexID] = intersect(Cortex_Incommon, Mice_STRUCT_T(mouseidx).REGION_NAME) ;

    CortexID_Mouse(mouseidx,:)    = CortexID;

end


[~, sortIdx]                      = sort(CortexID_Mouse(1, :) );

CortexID_Mouse_Sorted             = CortexID_Mouse(:, sortIdx );

Cortex_Incommon_Sorted            = Cortex_Incommon(sortIdx   );
 

%--------------------------------------------------------------------------
REGION_NAME_T = {};

for mouseidx  = 1:5

    REGION_NAME_T = [REGION_NAME_T ; Mice_STRUCT_T(mouseidx).REGION_NAME] ;

end

REGION_NAME_T = unique(REGION_NAME_T);



path                            = 'D:\XM\FC_M94\';
addpath([path 'othercolor'],'begin')
Mycolor                         = othercolor('Cat_12');

color_inuse = [
    
    220,20,60;
    255,99,71;

    154,205,50;
    124,252,0;
    
    0,255,255;
    0,206,209;

    127,255,212;
    64,224,208;

    220,20,60;
    255,99,71;
    
    102,205,170;
    32,178,170; 

    250,128,114;
    30,144,255;

    233,150,122;
    255,215,0;

    138,43,226;
    123,104,238;

    100,149,237;
    0,191,255;

    32,178,170;
    0,191,255;
 
    188,143,143;
    147,112,219;

    123,104,238;
    255,105,180;
    255,20,147;
  
    ]/255;


REGION_COLOR_T    = zeros( size(REGION_NAME_T,1) , 3 );
REGION_COLOR_mark = zeros( size(REGION_NAME_T,1) , 1 );

for regionidx     = 1:size( Cortex_Incommon_Sorted,1 )

    regionindex   = find(  strcmp(  Cortex_Incommon_Sorted{regionidx},REGION_NAME_T  )  );
    REGION_COLOR_mark(regionindex,:) = 1;
    REGION_COLOR_T(regionindex,:)    = color_inuse(regionidx,:);

end

REGION_COLOR_T(find(1-REGION_COLOR_mark) , :) = Mycolor(find(1-REGION_COLOR_mark)*1,:);

for mouseidx=1:5

    for regionidx   = 1:size( Mice_STRUCT_T(mouseidx).REGION_NAME,1)

    Mice_STRUCT_T(mouseidx).COLOR(regionidx,:) = REGION_COLOR_T(   find(  strcmp( Mice_STRUCT_T(mouseidx).REGION_NAME{regionidx} , REGION_NAME_T )  )  , :  );
  
    end

  mouseidx

end



figure(200)

for mouseidx                       = 1:5

CueShock_marker_color(Mice_STRUCT_T(5).NEURON.CUE_Frame(1,:) , Mice_STRUCT_T(5).NEURON.CUE_Frame(2,:) , Mice_STRUCT_T(5).NEURON.SHOCK_Frame(1,:) , Mice_STRUCT_T(5).NEURON.SHOCK_Frame(2,:) , 27 );

CortexID_Mouse_Sorted(mouseidx,:)

Regions_GROUP_Index

[~, index_2]                   = ismember(Regions_GROUP_Index, CortexID_Mouse_Sorted(mouseidx,:));

index_2_Regroup = nonzeros(index_2)

Trace_Region_Avg = []

for regionidx                  = 1:size(CortexID_Mouse_Sorted,2)

    Trace_Region_Avg(regionidx,:,mouseidx) = Mice_STRUCT_T(mouseidx).NEURON.Trace_unify_Region_Avg(Regions_GROUP_Index(regionidx) , : );

    Trace_Region_Avg_Norm      = rescale(Trace_Region_Avg( regionidx,:,mouseidx));

%     Trace_Region_Avg_Norm_T( regionidx,:,mouseidx) = Trace_Region_Avg_Norm;

%     plot(  Trace_Region_Avg_Norm + 27-(regionidx)+1 , 'LineWidth',0.6,'Color', color_inuse(index_2_Regroup(regionidx),:)  )
% 
%     hold on;

end
hold off 
end


%-----------------------------------------
%-----------------------------------------

Mouse_chosen=1



figure(201)
imshow(  uint8( Mice_STRUCT_T(Mouse_chosen).ATLAS2COMET )  );

figure(202)


imshow(uint8(255* ones(size(Mice_STRUCT_T(Mouse_chosen).ATLAS2COMET,1),size(Mice_STRUCT_T(Mouse_chosen).ATLAS2COMET,2)) ));
hold on;


for regionidx=1:length(Mice_STRUCT_T(Mouse_chosen).REGION_ID)

    BW = imbinarize(double(Mice_STRUCT_T(Mouse_chosen).ATLAS2COMET ==regionidx));
    boundaries = bwboundaries(BW);


    if (~isempty(boundaries))
        for k = 1:1

            b = boundaries{k};
            plot(b(:,2), b(:,1), 'Color' , [0 0 0], 'LineWidth', 1); % Plot boundary with green line
            %     fill(b(:,2)+1, b(:,1)+1, Mice_STRUCT_T(1).COLOR(11,:),'FaceAlpha',0.5,'LineStyle',':');
            hold on

        end
        hold on
    end
    hold on
    regionidx
end


hold on



%-----chosen Mouse=M94, for images--------------------------------
regionmoter = [9 10 31 32 ];%9 10 13 14 15 16 31 32 33 34 61 62 65 66 69 70 75 76]%
%-----chosen Mouse=5,M94, for images--------------------------------
for regionmoteridx=1:length(regionmoter)

BW         = imbinarize(double(Mice_STRUCT_T(Mouse_chosen).ATLAS2COMET ==regionmoter(regionmoteridx)));
boundaries = bwboundaries(BW);

for k = 1:1

    b = boundaries{k};
  
        fill(b(:,2)+1, b(:,1)+1, Mice_STRUCT_T(Mouse_chosen).COLOR(regionmoter(regionmoteridx),:),'FaceAlpha',0.5,'LineStyle',':');
    hold on

end

hold on

end



figure(204)
imshow(uint8(255* ones(size(Mice_STRUCT_T(Mouse_chosen).ATLAS2COMET,1),size(Mice_STRUCT_T(Mouse_chosen).ATLAS2COMET,2)) ));
hold on;
for regionidx=1:length(Mice_STRUCT_T(Mouse_chosen).REGION_ID)

    BW = imbinarize(double(Mice_STRUCT_T(Mouse_chosen).ATLAS2COMET ==regionidx));
    boundaries = bwboundaries(BW);


    if (~isempty(boundaries))
        for k = 1:1

            b = boundaries{k};
            plot(b(:,2), b(:,1), 'Color' , [0 0 0], 'LineWidth', 1); % Plot boundary with green line
            %     fill(b(:,2)+1, b(:,1)+1, Mice_STRUCT_T(1).COLOR(11,:),'FaceAlpha',0.5,'LineStyle',':');
            hold on

        end
        hold on
    end
    hold on
    regionidx
end

hold on

%-----chosen Mouse=M94, for images--------------------------------
regionmoter = [7 8 29 30 43 44]%[31 32 45 46]
%-----chosen Mouse=M94, for images--------------------------------

for regionmoteridx=1:length(regionmoter)

BW = imbinarize(double(Mice_STRUCT_T(Mouse_chosen).ATLAS2COMET ==regionmoter(regionmoteridx)));
boundaries = bwboundaries(BW);

for k = 1:1

    b = boundaries{k};
  
        fill(b(:,2)+1, b(:,1)+1, Mice_STRUCT_T(Mouse_chosen).COLOR(regionmoter(regionmoteridx),:),'FaceAlpha',0.5,'LineStyle',':');
    hold on

end

hold on

end

%%=========================================================================

close all

for mouseidx=1:5

    FIG_NEURON = figure(mouseidx)
    imshow(uint8(255* ones(size(Mice_STRUCT_T(mouseidx).ATLAS2COMET,1),size(Mice_STRUCT_T(mouseidx).ATLAS2COMET,2)) ));
    hold on;
    for regionidx=1:length(Mice_STRUCT_T(mouseidx).REGION_ID)

        BW = imbinarize(double(Mice_STRUCT_T(mouseidx).ATLAS2COMET ==regionidx));
        boundaries = bwboundaries(BW);

        if (~isempty(boundaries))
            for k = 1:1

                b = boundaries{k};
                plot(b(:,2), b(:,1), 'Color' , [0 0 0], 'LineWidth', 1); % Plot boundary with green line
                %     fill(b(:,2)+1, b(:,1)+1, Mice_STRUCT_T(1).COLOR(11,:),'FaceAlpha',0.5,'LineStyle',':');
                hold on

            end
            hold on
        end
        hold on
        regionidx
    end

    hold on
    scatter( Mice_STRUCT_T(mouseidx).COORDINATES(:,1) , Mice_STRUCT_T(mouseidx).COORDINATES(:,2),8,Mice_STRUCT_T(mouseidx).COLOR(Mice_STRUCT_T(mouseidx).Neuron_REGION_ID,:) ,'filled' , 'MarkerEdgeColor', [0.5 0.5 0.5], 'LineWidth', 0.1)
    hold off

    FIG_COLORBAR=figure(50+mouseidx);
    FIG_COLORBAR.Position = [500,50,180,900]
    COLOR_MOUSE = Mice_STRUCT_T(mouseidx).COLOR(unique(Mice_STRUCT_T(mouseidx).Neuron_REGION_ID) , :);
    scatter( repmat(1,1,size(COLOR_MOUSE,1)), 1:size(COLOR_MOUSE,1) , 100,COLOR_MOUSE ,'filled')

    text(repmat(1.15,1,size(COLOR_MOUSE,1)),[1:size(COLOR_MOUSE,1)]-0.1,Mice_STRUCT_T(mouseidx).REGION_NAME(unique(Mice_STRUCT_T(mouseidx).Neuron_REGION_ID)) , 'FontSize', 12 )
    axis off;

    %pause

end



%---------------Chosen certain cortical regions----------------------------
%-----chosen Mouse=M94, for images--------------------------------
% regionmoter = [9 10 13 14 15 16 31 32 33 34 61 62 65 66 69 70 75 76]%
regionmoter = [9 10 31 32 ]%
%-----chosen Mouse=5,M94, for images--------------------------------
for mouseidx=1

    indicesInA = find(ismember(Mice_STRUCT_T(mouseidx).Neuron_REGION_ID, regionmoter));

    REGION_ID_chosen  = Mice_STRUCT_T(mouseidx).Neuron_REGION_ID;

    FIG_NEURON = figure(60+mouseidx)
    imshow(uint8(255* ones(size(Mice_STRUCT_T(mouseidx).ATLAS2COMET,1),size(Mice_STRUCT_T(mouseidx).ATLAS2COMET,2)) ));
    hold on;
    for regionidx=1:length(Mice_STRUCT_T(mouseidx).REGION_ID)

        BW = imbinarize(double(Mice_STRUCT_T(mouseidx).ATLAS2COMET ==regionidx));
        boundaries = bwboundaries(BW);


        if (~isempty(boundaries))
            for k = 1:1

                b = boundaries{k};
                plot(b(:,2), b(:,1), 'Color' , [0 0 0], 'LineWidth', 1); % Plot boundary with green line
                %     fill(b(:,2)+1, b(:,1)+1, Mice_STRUCT_T(1).COLOR(11,:),'FaceAlpha',0.5,'LineStyle',':');
                hold on

            end
            hold on
        end
        hold on
        regionidx
    end

    hold on
    scatter( Mice_STRUCT_T(mouseidx).COORDINATES(indicesInA,1) , Mice_STRUCT_T(mouseidx).COORDINATES(indicesInA,2),12,Mice_STRUCT_T(mouseidx).COLOR(REGION_ID_chosen(indicesInA),:) ,'filled', 'MarkerEdgeColor', [0.5 0.5 0.5], 'LineWidth', 0.1 )
    hold off

    %pause

end



regionvisal = [7 8 29 30 43 44]%[31 32 45 46]

for mouseidx=1

    indicesInA = find(ismember(Mice_STRUCT_T(mouseidx).Neuron_REGION_ID, regionvisal));

    REGION_ID_chosen  = Mice_STRUCT_T(mouseidx).Neuron_REGION_ID;

    FIG_NEURON = figure(70+mouseidx)
    imshow(uint8(255* ones(size(Mice_STRUCT_T(mouseidx).ATLAS2COMET,1),size(Mice_STRUCT_T(mouseidx).ATLAS2COMET,2)) ));
    hold on;
    for regionidx=1:length(Mice_STRUCT_T(mouseidx).REGION_ID)

        BW = imbinarize(double(Mice_STRUCT_T(mouseidx).ATLAS2COMET ==regionidx));
        boundaries = bwboundaries(BW);

        if (~isempty(boundaries))
            for k = 1:1

                b = boundaries{k};
                plot(b(:,2), b(:,1), 'Color' , [0 0 0], 'LineWidth', 1); % Plot boundary with green line
                %     fill(b(:,2)+1, b(:,1)+1, Mice_STRUCT_T(1).COLOR(11,:),'FaceAlpha',0.5,'LineStyle',':');
                hold on

            end
            hold on
        end
        hold on
        regionidx
    end
    hold on
    scatter( Mice_STRUCT_T(mouseidx).COORDINATES(indicesInA,1) , Mice_STRUCT_T(mouseidx).COORDINATES(indicesInA,2),12,Mice_STRUCT_T(mouseidx).COLOR(REGION_ID_chosen(indicesInA),:) ,'filled', 'MarkerEdgeColor', [0.5 0.5 0.5], 'LineWidth', 0.1 )
    hold off

    %pause

end
%---------------Chosen certain cortical regions----------------------------
%%=========================================================================


%%=========================================================================
FIG_SHOCK     = figure(300)
plot( Mice_STRUCT_T(1).NEURON.SHOCK     )

FIG_SPON      = figure(301)
plot( Mice_STRUCT_T(1).NEURON.SPON      )

FIG_CUE       = figure(302)
plot( Mice_STRUCT_T(1).NEURON.CUE       )

FIG_CUE2SHOCK = figure(303)
plot( Mice_STRUCT_T(1).NEURON.CUE2SHOCK )

FIG_LEARN     = figure(304)
plot( Mice_STRUCT_T(1).NEURON.LEARN     )

%%-------------------------------------------------------------------------

Numframe                            = [];
NumFrame_min                        = 100000;
CUESHOCK_align                      = []
CUESHOCK_align_diff                 = []

for mouseidx                        = 1 : size(Mice_STRUCT_T,2)


    Numframe                        = [  Numframe size(Mice_STRUCT_T(mouseidx).NEURON.TIME , 2) ];
    SHOCK                           = Mice_STRUCT_T(mouseidx).NEURON.SHOCK;
    CUESHOCK_align(mouseidx,:)      = [Mice_STRUCT_T(mouseidx).NEURON.CUESHOCK_Frame , length(Mice_STRUCT_T(mouseidx).NEURON.TIME)];

    CUESHOCK_align_diff(mouseidx,:) = diff(CUESHOCK_align(mouseidx,:));
     
end

%%-------------------------------------------------------------------------

framemin               = min(CUESHOCK_align_diff,[],1);
framemindiff           = [2303 framemin]
framemindiff_cueshock  = framemindiff(2:end);
framemindiff_cueshock1 = ceil(framemindiff_cueshock*2/3)
framemindiff_cueshock2 = framemindiff_cueshock-framemindiff_cueshock1-1

Frame_Inuse_Mouse      = {}

for mouseidx           = 1 : size(Mice_STRUCT_T,2)


    CUESHOCK_chosen    = CUESHOCK_align(mouseidx,:);

    Frame_Inuse        = zeros(1,length(Mice_STRUCT_T(mouseidx).NEURON.TIME));

    Frame_Inuse(1,CUESHOCK_chosen(1) - framemindiff(1) + 1:CUESHOCK_chosen(1)-1) = 1;

    for trialidx       = 1:9

        Frame_Inuse(1,CUESHOCK_chosen(trialidx):CUESHOCK_chosen(trialidx)+framemindiff_cueshock1(trialidx)) = 1;
        Frame_Inuse(1,CUESHOCK_chosen(trialidx+1)-1:-1:CUESHOCK_chosen(trialidx+1)-framemindiff_cueshock2(trialidx)) = 1;

    end

    Frame_Inuse(1,CUESHOCK_chosen(9+1):CUESHOCK_chosen(9+1)+framemindiff_cueshock(9+1)) = 1;

    sum(Frame_Inuse)

    Frame_Inuse_Mouse{mouseidx}=Frame_Inuse;



end

%%-------------------------------------------------------------------------
%%-------------------------------------------------------------------------

%%==================Animation for activity of neurons during SHOCK=========
% % % % for mouseidx=1:5
% % % % 
% % % % 
% % % %     for trialidx=1:5
% % % % 
% % % %         ALPHA_T=[]
% % % %         Zeropad_left    = 10;
% % % %         Zeropad_right   = 10;
% % % % 
% % % %         SHOCK_frame     = find(Mice_STRUCT_T(mouseidx).NEURON.SHOCK(logical(Frame_Inuse_Mouse{mouseidx}))==trialidx);
% % % %         SHOCK_frame_z   = [SHOCK_frame(1)-Zeropad_left:SHOCK_frame(1)-1, SHOCK_frame, SHOCK_frame(end)+1:SHOCK_frame(end)+Zeropad_right];
% % % % 
% % % %         for frameidx    = 1:length(SHOCK_frame_z)
% % % % 
% % % %             FIG_NEURON  = figure(70+mouseidx)
% % % %             imshow(uint8(255*rot90(Mice_STRUCT_T(mouseidx).ATLAS90,1)))
% % % %             hold on
% % % % 
% % % %             Neuron_ID   = Mice_STRUCT_T(mouseidx).Neuron_REGION_ID   ;
% % % %             COLOR       = Mice_STRUCT_T(mouseidx).COLOR(Neuron_ID,:) ;
% % % %             %             ALPHA_T     = Mice_STRUCT_T(mouseidx).NEURON.Trace_unify(:,logical(Frame_Inuse_Mouse{mouseidx}));
% % % %             %             ALPHA_T     = (ALPHA_T(:,SHOCK_frame_z(frameidx)));
% % % % 
% % % %             SHOCK_CHOSEN = find(Mice_STRUCT_T(1).NEURON.SHOCK==trialidx)
% % % %             ALPHA_T     = Mice_STRUCT_T(mouseidx).NEURON.Trace_unify(: , SHOCK_CHOSEN(1)-10:SHOCK_CHOSEN(end)+10 );
% % % %             ALPHA = ALPHA_T(:,frameidx);
% % % % 
% % % % 
% % % %             %             COLOR_TRACE = COLOR .* ;
% % % % 
% % % %             h=scatter( Mice_STRUCT_T(mouseidx).COORDINATES(:,1) , Mice_STRUCT_T(mouseidx).COORDINATES(:,2) , 8 , COLOR , 'filled' );
% % % %             h.MarkerFaceAlpha = 'flat'; % 设置标记的面透明度为可变
% % % %             h.AlphaData = ALPHA;  % 为每个点设置透明度
% % % %             hold off
% % % % 
% % % %             FRAME       = getframe;
% % % % 
% % % %             imwrite(  FRAME.cdata , [Dir  Mice_STRUCT_T(mouseidx).NAME '_SHOCK-TRIAL_' num2str(trialidx) '1.tiff'] , "WriteMode" , "append"  );
% % % % 
% % % %         end
% % % % 
% % % %         %         %pause
% % % % 
% % % %     end
% % % % 
% % % %     %pause
% % % % 
% % % % 
% % % % end

%%-------------------------------------------------------------------------
% % % %    
% % % % Neuron_ID           = Mice_STRUCT_T(mouseidx).Neuron_REGION_ID   ;
% % % %     COLOR               = Mice_STRUCT_T(mouseidx).COLOR(Neuron_ID,:) ;
% % % % 
% % % %        ALPHA_T          = []
% % % %         Zeropad_left    = 10;
% % % %         Zeropad_right   = 10;
% % % %     ALPHA_T             = Mice_STRUCT_T(mouseidx).NEURON.Trace_unify(:,logical(Frame_Inuse_Mouse{mouseidx}));
% % % %   
% % % % 
% % % % for frameidx            = 1:39
% % % %     
% % % %     ALPHA_T_MAX         = [];
% % % % 
% % % %     for trialidx        = 1:5
% % % % 
% % % %         SHOCK_frame     = find(Mice_STRUCT_T(mouseidx).NEURON.SHOCK(logical(Frame_Inuse_Mouse{mouseidx}))==trialidx);
% % % %         SHOCK_frame_z   = [SHOCK_frame(1)-Zeropad_left:SHOCK_frame(1)-1, SHOCK_frame, SHOCK_frame(end)+1:SHOCK_frame(end)+Zeropad_right];
% % % % 
% % % %         ALPHA_T         = ALPHA_T(:,SHOCK_frame_z(frameidx));
% % % %         ALPHA_T_MAX     = cat(2,ALPHA_T_MAX,ALPHA_T);
% % % % 
% % % %     end
% % % % 
% % % % 
% % % %     FIG_NEURON          = figure(80)
% % % %     imshow( uint8(255*rot90(Mice_STRUCT_T(mouseidx).ATLAS90,1)) )
% % % %     hold on
% % % % 
% % % %    
% % % %         h           = scatter( Mice_STRUCT_T(mouseidx).COORDINATES(:,1) , Mice_STRUCT_T(mouseidx).COORDINATES(:,2) , 8 , COLOR , 'filled' );
% % % %     h.MarkerFaceAlpha = 'flat'; % 设置标记的面透明度为可变
% % % %     h.AlphaData = ALPHA_T;  % 为每个点设置透明度
% % % %     hold off
% % % % 
% % % %     FRAME       = getframe;
% % % % 
% % % %     imwrite(  FRAME.cdata , [Dir  Mice_STRUCT_T(mouseidx).NAME '_SHOCK-TRIAL_' num2str(trialidx) '.tiff'] , "WriteMode" , "append"  );
% % % % 
% % % %     %         %pause
% % % % 
% % % % end
%%==================Animation for activity of neurons during SHOCK=========




%%-------------------------------------------------------------------------
%%-------------------------------------------------------------------------

figure(400)
    plot( Frame_Inuse_Mouse{1}/2)
    hold on;
    plot(Mice_STRUCT_T(1).NEURON.SHOCK)
    hold on;
    plot(Mice_STRUCT_T(1).NEURON.CUE)

figure(401)
    for mouseidx=1:5
        plot(Mice_STRUCT_T(mouseidx).NEURON.CUE(logical(Frame_Inuse_Mouse{mouseidx})));
        hold on
        plot(Mice_STRUCT_T(mouseidx).NEURON.SHOCK(logical(Frame_Inuse_Mouse{mouseidx})));
    end

%%-------------------------------------------------------------------------

mouseidx                 = 1
miniscopeshockbeginFrame = find(diff(Mice_STRUCT_T(mouseidx).NEURON.SHOCK(logical(Frame_Inuse_Mouse{mouseidx})))>0)+1
miniscopeshockendFrame   = find(diff(Mice_STRUCT_T(mouseidx).NEURON.SHOCK(logical(Frame_Inuse_Mouse{mouseidx})))<0)
miniscopecuebeginFrame   = find(diff(Mice_STRUCT_T(mouseidx).NEURON.CUE(logical(Frame_Inuse_Mouse{mouseidx})))>0)+1
miniscopecueendFrame     = find(diff(Mice_STRUCT_T(mouseidx).NEURON.CUE(logical(Frame_Inuse_Mouse{mouseidx})))<0)

%%-------------------------------------------------------------------------

close all

Trace_Region_Avg        = []


for mouseidx                    = 1:5;%1:5


    Trials                      = Mice_STRUCT_T(mouseidx).NEURON.SHOCK(logical(Frame_Inuse_Mouse{mouseidx})) + Mice_STRUCT_T(mouseidx).NEURON.CUE(logical(Frame_Inuse_Mouse{mouseidx}) );

    Trials_diff                 = [(diff(Trials)),0];

    result(mouseidx,:)          = cumtrapz(1:length(Trials), Trials_diff>0)+1;

    %%--------------------------------------
%     FIG_TRIAL = figure(500)
%     plot( Trials,'k')
%     hold on;
%     plot( result(mouseidx,:),'b')
    %%--------------------------------------

% 
%     for regionidx=1:size(Mice_STRUCT_T(mouseidx).REGION_ID,1)
% 
%         Mice_STRUCT_T(mouseidx).NEURON.Trace_unify   =  rescale(Mice_STRUCT_T(mouseidx).NEURON.Trace')';
% 
%         Mice_STRUCT_T(mouseidx).NEURON.Trace_unify_Region_Avg(regionidx,:) = mean( Mice_STRUCT_T(mouseidx).NEURON.Trace_unify( Mice_STRUCT_T(mouseidx).Neuron_REGION_ID==regionidx , : ) , 1 );
% 
%     end


    Trace_unify_SUM(mouseidx,:)  = sum( Mice_STRUCT_T(mouseidx).NEURON.Trace_unify(:,logical(Frame_Inuse_Mouse{mouseidx})) , 1 );
    Trace_unify_NUM(mouseidx,:)  = size( Mice_STRUCT_T(mouseidx).NEURON.Trace_unify , 1 );
    Trace_unify_AVG(mouseidx,:)  = mean( Mice_STRUCT_T(mouseidx).NEURON.Trace_unify(:,logical(Frame_Inuse_Mouse{mouseidx})) , 1 );

    Cortex_unify_AVG(mouseidx,:) = mean( Mice_STRUCT_T(mouseidx).REGION_ACT(:,logical(Frame_Inuse_Mouse{mouseidx})) , 1 );

    RawAct_unify_AVG(mouseidx,:) = mean( Mice_STRUCT_T(mouseidx).NEURON.RawAct_unify(:,logical(Frame_Inuse_Mouse{mouseidx})) , 1 );

    SPEED_unify(mouseidx,:)      = Mice_STRUCT_T(mouseidx).NEURON.Speed(logical(Frame_Inuse_Mouse{mouseidx}))/max(Mice_STRUCT_T(mouseidx).NEURON.Speed(logical(Frame_Inuse_Mouse{mouseidx}))) ;
    SPEED_unify(mouseidx,~isnan(SPEED_unify(mouseidx,:))) = mapminmax( SPEED_unify(mouseidx,~isnan(SPEED_unify(mouseidx,:))) , 0 ,1  );

    %%--------------------------------------
%     FIG_CORTEX_TRACE = figure(501)
% 
%     plot(  Cortex_unify_AVG(mouseidx , : )   , 'LineWidth',1,'Color', 'k'  )
%     hold on;
%     plot(  3*Trace_unify_AVG( mouseidx , : ) , 'LineWidth',1,'Color', 'r'  )
%     hold on;
%     plot(  SPEED_unify( mouseidx , : )       , 'LineWidth',1,'Color', 'b'  )
    %%--------------------------------------

%     FIG_REGION_TRACE = figure(502)
%     %  Hobj = gca;
%     %  set(Hobj,'position',[100,100,200,1000])
%     CueShock_marker_color(miniscopecuebeginFrame, miniscopecueendFrame,miniscopeshockbeginFrame,miniscopeshockendFrame,27);





    for regionidx                = 1:size(CortexID_Mouse_Sorted,2)


        Trace_Region_Avg(regionidx,:,mouseidx)  = Mice_STRUCT_T(mouseidx).NEURON.Trace_unify_Region_Avg(CortexID_Mouse_Sorted(mouseidx,regionidx) , logical(Frame_Inuse_Mouse{mouseidx}));
        

        RawAct_Region_Avg(regionidx,:,mouseidx)  = Mice_STRUCT_T(mouseidx).NEURON.RawAct_uniTraceAVG(CortexID_Mouse_Sorted(mouseidx,regionidx) , logical(Frame_Inuse_Mouse{mouseidx}));
        

        Neuron_Region_Num(regionidx,:,mouseidx) = Mice_STRUCT_T(mouseidx).REGION_ID( CortexID_Mouse_Sorted(mouseidx,regionidx) );

        Cortex_Region_Avg(regionidx,:,mouseidx) = Mice_STRUCT_T(mouseidx).REGION_ACT(CortexID_Mouse_Sorted(mouseidx,regionidx),logical(Frame_Inuse_Mouse{mouseidx}));


%         FIG_REGION_TRACE = figure(502)
%         Hobj = gca;

        %         plot(  Cortex_Region_Avg(regionidx,:,mouseidx) + (regionidx)-1 , 'LineWidth',1,'Color', 'k'  )
        %         hold on;
        Trace_Region_Avg_Norm = rescale(Trace_Region_Avg( regionidx,:,mouseidx));
        Trace_unify_AVG_Norm = rescale(Trace_unify_AVG( mouseidx,:));

%         plot(  Trace_Region_Avg_Norm + 27-(regionidx)+1 , 'LineWidth',2,'Color', color_inuse(regionidx,:)  )
%         hold on;
%         plot(Trace_unify_AVG_Norm+ 27-(regionidx)+1 , 'LineWidth',2,'Color', 'k'  )
% 
%         hold on;



        %         plot(  (Mice_STRUCT_T(mouseidx).NEURON.SHOCK(logical(Frame_Inuse_Mouse{mouseidx}))>0)+ (regionidx)-1,'k'  )
        %         hold on;
        %         plot(  (Mice_STRUCT_T(mouseidx).NEURON.CUE(logical(Frame_Inuse_Mouse{mouseidx}))>0)+ (regionidx)-1,'k'    )
        %         hold on;
%         text(7000,27-regionidx+1, Cortex_Incommon_Sorted{regionidx});
%         set(Hobj.XAxis,'Visible','off')
%         set(gca, 'Color', 'none');            % 轴背景透明
% 
%         %         set(gcf, 'Color', 'none');            % 图形背景透明
%         grid off
%         box off
% 
%         %%--------------------------------------
%         FIG_CORTEX_ALL = figure(503)
%         title(Cortex_Incommon_Sorted{regionidx})
%         plot(  Cortex_Region_Avg(regionidx , : , mouseidx)' + (regionidx)-1 ,'Color', color_inuse(regionidx,:)  )
%         hold on
%         plot(  Cortex_unify_AVG(mouseidx , : ) + (regionidx)-1  , 'LineWidth',1,'Color', 'k'  )
%         hold on
%         plot(  (Mice_STRUCT_T(mouseidx).NEURON.SHOCK(logical(Frame_Inuse_Mouse{mouseidx}))>0)+ (regionidx)-1,'k'  )
%         hold on;
%         plot(  (Mice_STRUCT_T(mouseidx).NEURON.CUE(logical(Frame_Inuse_Mouse{mouseidx}))>0)+ (regionidx)-1,'k'    )
%         hold on;
% 
%         text(6900,1*regionidx-1, Cortex_Incommon_Sorted{regionidx});
%         %         %pause
%         %%--------------------------------------
% 
% 
%         %%--------------------------------------
%         FIG_CORTEX_ALL = figure(504)
%         title(Cortex_Incommon_Sorted{regionidx})
% 
%         plot(  (Cortex_Region_Avg(regionidx , : , mouseidx)')  ,'Color', color_inuse(regionidx,:)  )
%         hold on
%         plot(  (Cortex_unify_AVG(mouseidx , : ))   , 'LineWidth',1,'Color', 'k'  )
% 
% 
%         %         %pause
%         %%--------------------------------------
% 
%         %%--------------------------------------
%         FIG_TRACE_ALL = figure(505+mouseidx)
%         title(Cortex_Incommon_Sorted{regionidx})
% 
%         plot(  3*(Trace_Region_Avg(regionidx , : , mouseidx)')+ (regionidx)-1 , 'LineWidth',2 ,'Color', color_inuse(regionidx,:)  )
%         hold on
%         plot(  3*(Trace_unify_AVG(mouseidx , : ))+ (regionidx)-1   , 'LineWidth',2,'Color', 'k'  )
%         hold on
%         plot(  (Mice_STRUCT_T(mouseidx).NEURON.SHOCK(logical(Frame_Inuse_Mouse{mouseidx}))>0)+ (regionidx)-1,'k'  )
%         hold on;
%         plot(  (Mice_STRUCT_T(mouseidx).NEURON.CUE(logical(Frame_Inuse_Mouse{mouseidx}))>0)+ (regionidx)-1,'k'    )
%         hold on;
% 
%         text(6900,1*regionidx-1, Cortex_Incommon_Sorted{regionidx});
% 

        %         %pause
        %%--------------------------------------


    end




    disp(   ['Mouse ID: ' num2str(mouseidx)]   );

    %     %pause


end




for mouseidx=1:5
figure(601)

subplot(5,1,mouseidx)
CueShock_marker_color_Y(miniscopecuebeginFrame, miniscopecueendFrame,miniscopeshockbeginFrame,miniscopeshockendFrame,0.15);

% for mouseidx=1:5
%   plot(  (Trace_unify_AVG(mouseidx,:))   , 'LineWidth',1  )
% hold on  
% end

plot(  Trace_unify_AVG(mouseidx,:)   , 'LineWidth',2,'Color', 'k'  )
hold on
Trace_Region_Avg_Tmouse = Trace_Region_Avg(:,:,mouseidx);
Trace_Region_Avg_Tmouse_MO =  Trace_Region_Avg_Tmouse([3,4,11,12],:);
plot(  mean(Trace_Region_Avg_Tmouse_MO,1)   , 'LineWidth',2,'Color', 'r'  )

figure(602)
subplot(5,1,mouseidx)
CueShock_marker_color_Y(miniscopecuebeginFrame, miniscopecueendFrame,miniscopeshockbeginFrame,miniscopeshockendFrame,0.15);
plot(  Trace_unify_AVG(mouseidx,:)   , 'LineWidth',2,'Color', 'k'  )
hold on
Trace_Region_Avg_Tmouse = Trace_Region_Avg(:,:,mouseidx);
Trace_Region_Avg_Tmouse_RSP = Trace_Region_Avg_Tmouse([5,6],:) ;
plot(  mean(Trace_Region_Avg_Tmouse_RSP,1)   , 'LineWidth',2,'Color', 'r'  )


figure(603)
subplot(5,1,mouseidx)
CueShock_marker_color_Y(miniscopecuebeginFrame, miniscopecueendFrame,miniscopeshockbeginFrame,miniscopeshockendFrame,0.15);
plot(  Trace_unify_AVG(mouseidx,:)   , 'LineWidth',2,'Color', 'k'  )
hold on
Trace_Region_Avg_Tmouse = Trace_Region_Avg(:,:,mouseidx);
Trace_Region_Avg_Tmouse_SSP = Trace_Region_Avg_Tmouse([7,8,17:24],:);
plot(  mean(Trace_Region_Avg_Tmouse_SSP,1)   , 'LineWidth',2,'Color', 'r'  )



figure(604)
subplot(5,1,mouseidx)
CueShock_marker_color_Y(miniscopecuebeginFrame, miniscopecueendFrame,miniscopeshockbeginFrame,miniscopeshockendFrame,0.15);
plot(  Trace_unify_AVG(mouseidx,:)   , 'LineWidth',2,'Color', 'k'  )
hold on
Trace_Region_Avg_Tmouse = Trace_Region_Avg(:,:,mouseidx);
Trace_Region_Avg_Tmouse_VIS = Trace_Region_Avg_Tmouse([1,2,9,10,15,16],:) ;
plot(  mean(Trace_Region_Avg_Tmouse_VIS,1)   , 'LineWidth',2,'Color', 'r'  )


end


for regionidx = 1:27

    handle = figure(700)

    for mouseidx=1:5

        subplot(5,1,mouseidx)
        CueShock_marker_color_Y(miniscopecuebeginFrame, miniscopecueendFrame,miniscopeshockbeginFrame,miniscopeshockendFrame,0.15);

        plot(  Trace_unify_AVG(mouseidx,:)   , 'LineWidth',2,'Color', 'k'  )
        hold on
        Trace_Region_Avg_Tmouse = Trace_Region_Avg(:,:,mouseidx);
        Trace_Region_Avg_Tmouse_chosen =  Trace_Region_Avg_Tmouse(regionidx,:);
        plot(  mean(Trace_Region_Avg_Tmouse_chosen,1)   , 'LineWidth',2,'Color', 'r'  )
        t1 = title(Regions_overlap_inseq{regionidx})

    end

%     %pause

    close(handle)

end



 Trace_Region_Avg_TmouseAvg = mean( Trace_Region_Avg , 3 );

close all



 for regionidx = 1:13%[3 6 7 8 9 13 12]%


     handle          = figure(700+regionidx)
     handle.Position = [100 100 1000 200];

     CueShock_marker_color_Y(miniscopecuebeginFrame, miniscopecueendFrame,miniscopeshockbeginFrame,miniscopeshockendFrame,0.15);

     plot(  mean( Trace_Region_Avg_TmouseAvg , 1)   , 'LineWidth',1,'Color', 'k'  )
     hold on

     Trace_Region_Avg_TmouseAvg_chosen = Trace_Region_Avg_TmouseAvg([2*regionidx-1 2*regionidx],:);
     plot(  mean(Trace_Region_Avg_TmouseAvg_chosen,1)   , 'LineWidth',1,'Color', 'r'  )
     hold on

     temp = mean( Trace_Region_Avg_TmouseAvg , 1);
     x = [0, miniscopecuebeginFrame(2)]; % x 轴起点和终点
     y = [min(temp(1:miniscopecuebeginFrame(2))), min(temp(1:miniscopecuebeginFrame(2)))]; % 水平线的 y 值
     line(x, y, 'Color', [0/255, 0/255, 204/255], 'LineStyle', ':', 'LineWidth', 2); % 绘制绿色点线
     hold on

     x = [miniscopecuebeginFrame(2) length(temp)]; % x 轴起点和终点
     y = [min(temp(miniscopecuebeginFrame(2):end)), min(temp(miniscopecuebeginFrame(2):end))]; % 水平线的 y 值
     line(x, y, 'Color', [0/255, 0/255, 204/255], 'LineStyle', ':', 'LineWidth', 2); % 绘制绿色点线
     hold on

     ylim([0.04 0.14])
     set(gca, 'XColor', 'none');
     set(gca, 'XTickLabel', []);
     set(gca, 'YColor', 'none');
     set(gca, 'YTickLabel', []);
     t1 = title(Regions_overlap_inseq{2*regionidx-1})
     

%      %pause

%      close(handle)

 end


     handle = figure(777)
     handle.Position = [100 100 800 200];

     CueShock_marker_color_Y(miniscopecuebeginFrame, miniscopecueendFrame,miniscopeshockbeginFrame,miniscopeshockendFrame,0.15);

     plot(  mean( Trace_Region_Avg_TmouseAvg , 1)   , 'LineWidth',1,'Color', 'k'  )
     hold on

   
     temp = mean( Trace_Region_Avg_TmouseAvg , 1);
     x = [0, miniscopecuebeginFrame(2)]; % x 轴起点和终点
     y = [min(temp(1:miniscopecuebeginFrame(2))), min(temp(1:miniscopecuebeginFrame(2)))]; % 水平线的 y 值
     line(x, y, 'Color', [0/255, 0/255, 204/255], 'LineStyle', ':', 'LineWidth', 2); % 绘制绿色点线
     hold on

     x = [miniscopecuebeginFrame(2) length(temp)]; % x 轴起点和终点
     y = [min(temp(miniscopecuebeginFrame(2):end)), min(temp(miniscopecuebeginFrame(2):end))]; % 水平线的 y 值
     line(x, y, 'Color', [0/255, 0/255, 204/255], 'LineStyle', ':', 'LineWidth', 2); % 绘制绿色点线
     hold on

     ylim([0.04 0.14])
     set(gca, 'XColor', 'none');
     set(gca, 'XTickLabel', []);
     set(gca, 'YColor', 'none');
     set(gca, 'YTickLabel', []);
     t1 = title(Regions_overlap_inseq{2*regionidx-1})
     

%      %pause

%      close(handle)

%  end


Trace_Region_Avg_TmouseAvg_avg = mean( Trace_Region_Avg_TmouseAvg , 1);

for trialidx=1:5

%     for mouseidx=1:5
%         SHOCK_frame   = find(Mice_STRUCT_T(mouseidx).NEURON.SHOCK(logical(Frame_Inuse_Mouse{mouseidx}))==trialidx);
% 
%         SHOCK_frame_z = [SHOCK_frame(1)-Zeropad_left:SHOCK_frame(1)-1, SHOCK_frame, SHOCK_frame(end)+1:SHOCK_frame(end)+Zeropad_right];

        Trace_Region_Avg_TmouseAvg_avg_trial(:,trialidx) = Trace_Region_Avg_TmouseAvg_avg(miniscopeshockbeginFrame(trialidx)-5:miniscopeshockbeginFrame(trialidx)+69);
%     end


end


Trace_Region_Avg_TmouseAvg_avg_trial_avg = mean(Trace_Region_Avg_TmouseAvg_avg_trial,2);


close all


 for regionidx = [3 6 7 8 9]%1:13

     handle          = figure(720+regionidx)
     handle.Position = [100+regionidx*100 100 100 200];

       t1 = 6; % 开始时间
        t2 = 6+20; % 结束时间
        % 添加阴影
        px = [t1, t2, t2, t1]; % 阴影的 x 坐标
        py = [0, 0, 1, 1]; % 阴影的 y 坐标
        patch(px, py, [178 102 255]/255, 'FaceAlpha', 0.6, 'EdgeColor', 'none'); % 添加阴影
        hold on;

     plot(  Trace_Region_Avg_TmouseAvg_avg_trial_avg   , 'LineWidth',2,'Color', 'k'  )
     hold on

     Trace_Region_Avg_TmouseAvg_chosen = mean( Trace_Region_Avg_TmouseAvg([2*regionidx-1 2*regionidx],:) , 1);

    for trialidx=1:5

    Trace_Region_Avg_TmouseAvg_chosen_trial(:,trialidx) = Trace_Region_Avg_TmouseAvg_chosen(miniscopeshockbeginFrame(trialidx)-5:miniscopeshockbeginFrame(trialidx)+69);

end


Trace_Region_Avg_TmouseAvg_chosen_trial_avg = mean(Trace_Region_Avg_TmouseAvg_chosen_trial,2);

     
     plot(  (Trace_Region_Avg_TmouseAvg_chosen_trial_avg)   , 'LineWidth',2,'Color', 'r'  )
     hold on

     ylim([0.04 0.14])
     set(gca, 'XColor', 'none');
     set(gca, 'XTickLabel', []);
     set(gca, 'YColor', 'none');
     set(gca, 'YTickLabel', []);

%      t1            = title(Regions_overlap_inseq{2*regionidx-1})
   

 end


for mouseidx=1:5

   
    Timediff                      = [0;diff(Mice_STRUCT_T(mouseidx).NEURON.TIME)'];
    T_spk_diff(mouseidx,:)        = Timediff( logical(Frame_Inuse_Mouse{mouseidx}));

 
end




Frame_Learn_index = zeros(1 , length(Trace_unify_AVG_Norm));
Frame_Learn_index(1,1:miniscopecuebeginFrame(1)) = 1;
Frame_Learn_index(1,miniscopeshockendFrame(1):miniscopecuebeginFrame(2)) = 2;
Frame_Learn_index(1,miniscopeshockendFrame(2):miniscopecuebeginFrame(3)) = 3;
Frame_Learn_index(1,miniscopeshockendFrame(3):miniscopecuebeginFrame(4)) = 4;
Frame_Learn_index(1,miniscopeshockendFrame(4):miniscopecuebeginFrame(5)) = 5;
Frame_Learn_index(1,miniscopeshockendFrame(5):end) = 6;


diff_abs_learn = []

for mouseidx   = 1:5

Trace_Region_Avg_Tmouse = Trace_Region_Avg(:,:,mouseidx) ;

 for regionidx = 1:13


     handle = figure(720+regionidx)
     handle.Position = [100 100 800 200];

%      CueShock_marker_color_Y(miniscopecuebeginFrame, miniscopecueendFrame,miniscopeshockbeginFrame,miniscopeshockendFrame,0.15);

     Trace_Region_Avg_TmouseAvg_chosen = Trace_Region_Avg_Tmouse([2*regionidx-1 2*regionidx],:);
     diff_abs = mean(Trace_Region_Avg_TmouseAvg_chosen,1)  - mean( Trace_Region_Avg_Tmouse , 1) ;

 for learnidx=1:6

     diff_abs_learn(learnidx,mouseidx,regionidx) = mean(diff_abs(Frame_Learn_index==learnidx));
     diff_sum_learn(learnidx,mouseidx,regionidx) = sum(diff_abs(Frame_Learn_index==learnidx));

     T_spk_diff_learn             = T_spk_diff(mouseidx,Frame_Learn_index==learnidx) ;   
     T_spk_diff_learn_sum         = sum(T_spk_diff_learn)/1000;
 
     FR_diff_learn(learnidx,mouseidx,regionidx)                = diff_sum_learn(learnidx,mouseidx,regionidx)/T_spk_diff_learn_sum;


     plot(  mean(diff_abs(Frame_Learn_index==learnidx))  , 'LineWidth',1,'Color', 'k'  )
     hold on
     
%      Correlation_regions_mice(learnidx,mouseidx,regionidx) = xcorr2(mean(Trace_Region_Avg_TmouseAvg_chosen(:,Frame_Learn_index==learnidx),1));


 end
    
     set(gca, 'XColor', 'none');
     set(gca, 'XTickLabel', []);
     set(gca, 'YColor', 'none');
     set(gca, 'YTickLabel', []);
     t1 = title(Regions_overlap_inseq{2*regionidx-1})
     

%      %pause

%      close(handle)

 end

end



Frame_SHOCK_index = zeros( 1 , length(Trace_unify_AVG_Norm) );

Frame_SHOCK_index(1,miniscopeshockbeginFrame(1):miniscopeshockbeginFrame(1)+19) = 1;
Frame_SHOCK_index(1,miniscopeshockbeginFrame(2):miniscopeshockbeginFrame(2)+19) = 2;
Frame_SHOCK_index(1,miniscopeshockbeginFrame(3):miniscopeshockbeginFrame(3)+19) = 3;
Frame_SHOCK_index(1,miniscopeshockbeginFrame(4):miniscopeshockbeginFrame(4)+19) = 4;
Frame_SHOCK_index(1,miniscopeshockbeginFrame(5):miniscopeshockbeginFrame(5)+19) = 5;                   





diff_abs_SHOCK = []

for mouseidx   = 1:5

    Trace_Region_Avg_Tmouse = Trace_Region_Avg(:,:,mouseidx) ;

    for regionidx = 1:13

        Trace_Region_Avg_TmouseAvg_chosen = Trace_Region_Avg_Tmouse([2*regionidx-1 2*regionidx],:);
        diff_abs = mean(Trace_Region_Avg_TmouseAvg_chosen,1)  - mean( Trace_Region_Avg_Tmouse , 1) ;

        for learnidx=1:5

            diff_abs_SHOCK(learnidx,mouseidx,regionidx) = mean(diff_abs(Frame_SHOCK_index==learnidx));
       
            diff_sum_SHOCK(learnidx,mouseidx,regionidx) = sum(diff_abs(Frame_SHOCK_index==learnidx));

     T_spk_diff_SHOCK             = T_spk_diff(mouseidx,Frame_SHOCK_index==learnidx) ;   
     T_spk_diff_SHOCK_sum         = sum(T_spk_diff_SHOCK)/1000;
 
     FR_diff_SHOCK(learnidx,mouseidx,regionidx)                = diff_sum_SHOCK(learnidx,mouseidx,regionidx)/T_spk_diff_SHOCK_sum;


        end

      

    end

end


Frame_CUE_index = zeros( 1 , length(Trace_unify_AVG_Norm) );

Frame_CUE_index(1,miniscopecuebeginFrame(1):miniscopecuebeginFrame(1)+200) = 1;
Frame_CUE_index(1,miniscopecuebeginFrame(2):miniscopecuebeginFrame(2)+200) = 2;
Frame_CUE_index(1,miniscopecuebeginFrame(3):miniscopecuebeginFrame(3)+200) = 3;
Frame_CUE_index(1,miniscopecuebeginFrame(4):miniscopecuebeginFrame(4)+200) = 4;
Frame_CUE_index(1,miniscopecuebeginFrame(5):miniscopecuebeginFrame(5)+200) = 5;                   





diff_abs_SHOCK = []

for mouseidx   = 1:5

    Trace_Region_Avg_Tmouse = Trace_Region_Avg(:,:,mouseidx) ;

    for regionidx = 1:13

        Trace_Region_Avg_TmouseAvg_chosen = Trace_Region_Avg_Tmouse([2*regionidx-1 2*regionidx],:);
        diff_abs = mean(Trace_Region_Avg_TmouseAvg_chosen,1)  - mean( Trace_Region_Avg_Tmouse , 1) ;

        for learnidx=1:5

            diff_abs_CUE(learnidx,mouseidx,regionidx) = mean(diff_abs(Frame_CUE_index==learnidx));

            diff_sum_CUE(learnidx,mouseidx,regionidx) = sum(diff_abs(Frame_CUE_index==learnidx));

            T_spk_diff_CUE             = T_spk_diff(mouseidx,Frame_CUE_index==learnidx) ;
            T_spk_diff_CUE_sum         = sum(T_spk_diff_CUE)/1000;

            FR_diff_CUE(learnidx,mouseidx,regionidx)                = diff_sum_CUE(learnidx,mouseidx,regionidx)/T_spk_diff_CUE_sum;


        end


    end

end


%%========================Network analysis=================================

Frame_CUE_index = zeros( 1 , length(Trace_unify_AVG_Norm) );

cuelength = 200

Frame_CUE_index(1,miniscopecuebeginFrame(1):miniscopecuebeginFrame(1)+cuelength) = 1;
Frame_CUE_index(1,miniscopecuebeginFrame(2):miniscopecuebeginFrame(2)+cuelength) = 2;
Frame_CUE_index(1,miniscopecuebeginFrame(3):miniscopecuebeginFrame(3)+cuelength) = 3;
Frame_CUE_index(1,miniscopecuebeginFrame(4):miniscopecuebeginFrame(4)+cuelength) = 4;
Frame_CUE_index(1,miniscopecuebeginFrame(5):miniscopecuebeginFrame(5)+cuelength) = 5;        

% Frame_Learn_index = zeros(1 , length(Trace_unify_AVG_Norm));
% Frame_Learn_index(1,1:miniscopecuebeginFrame(1)) = 1;
% Frame_Learn_index(1,miniscopeshockendFrame(1):miniscopeshockendFrame(1)+200) = 2;
% Frame_Learn_index(1,miniscopeshockendFrame(2):miniscopeshockendFrame(2)+200) = 3;
% Frame_Learn_index(1,miniscopeshockendFrame(3):miniscopeshockendFrame(3)+200) = 4;
% Frame_Learn_index(1,miniscopeshockendFrame(4):miniscopeshockendFrame(4)+200) = 5;
% Frame_Learn_index(1,miniscopeshockendFrame(5):miniscopeshockendFrame(5)+200) = 6;

Frame_Learn_index = zeros(1 , length(Trace_unify_AVG_Norm));
Frame_Learn_index(1,1:miniscopecuebeginFrame(1)) = 1;
Frame_Learn_index(1,miniscopeshockendFrame(1):miniscopecuebeginFrame(2)) = 2;
Frame_Learn_index(1,miniscopeshockendFrame(2):miniscopecuebeginFrame(3)) = 3;
Frame_Learn_index(1,miniscopeshockendFrame(3):miniscopecuebeginFrame(4)) = 4;
Frame_Learn_index(1,miniscopeshockendFrame(4):miniscopecuebeginFrame(5)) = 5;
Frame_Learn_index(1,miniscopeshockendFrame(5):end) = 6;

Frame_chosen_index = Frame_Learn_index;




for mouseidx   = 1:5

Traces_all=[]
  
    Traces_all  = Mice_STRUCT_T(mouseidx).NEURON.Trace_unify(: , logical(Frame_Inuse_Mouse{mouseidx}));
Degree_centrality_Cells=[]

    for learnidx=1:max(Frame_chosen_index)
        
        Traces_learn  = Traces_all(:,Frame_chosen_index==learnidx);
     
        Correlation_Traces_learn  = corr(Traces_learn');
        Correlation_Traces_learn(isnan(Correlation_Traces_learn))=0;
     
         Degree_centrality_Cells(:, learnidx)  = (sum(Correlation_Traces_learn, 2)-1)/(size(Correlation_Traces_learn,1)-1);

    
    end

    Degree_centrality_Cells_Mice{mouseidx} = Degree_centrality_Cells;
mouseidx

end

for mouseidx=1:5

   figure(2000)
    imshow(uint8(255* ones(size(Mice_STRUCT_T(mouseidx).ATLAS2COMET,1),size(Mice_STRUCT_T(mouseidx).ATLAS2COMET,2)) ));
    hold on;
    for regionidx=1:length(Mice_STRUCT_T(mouseidx).REGION_ID)

        BW = imbinarize(double(Mice_STRUCT_T(mouseidx).ATLAS2COMET ==regionidx));
        boundaries = bwboundaries(BW);

        if (~isempty(boundaries))
            for k = 1:1

                b = boundaries{k};
                plot(b(:,2), b(:,1), 'Color' , [0 0 0], 'LineWidth', 1); % Plot boundary with green line
                %     fill(b(:,2)+1, b(:,1)+1, Mice_STRUCT_T(1).COLOR(11,:),'FaceAlpha',0.5,'LineStyle',':');
                hold on

            end
            hold on
        end
        hold on
        regionidx
    end



mouseidx
%pause

end

Bregma_Position = [985 973;
   933 973;
   1002 973;
  928 973;
 985 973 ]

color_mice = [1 0 0;1 0 1; 0 1 0; 0 0 1; 0 1 1]


figure(2000)
for mouseidx=1:5

    %     imshow(uint8(255* ones(size(Mice_STRUCT_T(mouseidx).ATLAS2COMET,1),size(Mice_STRUCT_T(mouseidx).ATLAS2COMET,2)) ));
    %     hold on;
    for regionidx=1:length(Mice_STRUCT_T(mouseidx).REGION_ID)

        BW = imbinarize(double(Mice_STRUCT_T(mouseidx).ATLAS2COMET ==regionidx));
        boundaries = bwboundaries(BW);

        if (~isempty(boundaries))

            for k = 1:1

                b = boundaries{k};
                b(:,2) = b(:,2) + 985 - Bregma_Position(mouseidx,1);
                plot(b(:,2), b(:,1), 'Color' , color_mice(mouseidx,:), 'LineWidth', 1); % Plot boundary with green line
                %     fill(b(:,2)+1, b(:,1)+1, Mice_STRUCT_T(1).COLOR(11,:),'FaceAlpha',0.5,'LineStyle',':');
                hold on

            end
            hold on
        end
        hold on
        regionidx
    end

    hold on

    %     scatter( Mice_STRUCT_T(mouseidx).COORDINATES(:,1) + 985 - Bregma_Position(mouseidx,1), Mice_STRUCT_T(mouseidx).COORDINATES(:,2),8,[1 1 0 ],'filled' )
    %     hold off
    % hold on

end


% 选择 colormap
colormap_name = 'jet'; % 例如选择 'jet' colormap
cmap = colormap(colormap_name); % 获取 colormap 的颜色矩阵

close all

for learnidx=1:6
 figure(2000+learnidx)
  imshow(uint8(255* ones(size(Mice_STRUCT_T(1).ATLAS2COMET,1),size(Mice_STRUCT_T(1).ATLAS2COMET,2)) ));
    hold on;
for mouseidx=[1 2 4 5]

Degree_centrality_Cells_color=[]
Degree_centrality_Cells_norm = rescale( Degree_centrality_Cells_Mice{mouseidx});






% 将灰度值映射到 colormap 的颜色
indices = round(Degree_centrality_Cells_norm(:,learnidx) * (size(cmap, 1) - 1)) + 1; % 映射到 [1, N] 的索引
indices = max(min(indices, size(cmap, 1)), 1); % 确保索引在有效范围内
Degree_centrality_Cells_color(:,:,learnidx) = cmap(indices, :); % 获取对应的颜色


  hold on
   
    for regionidx=1:length(Mice_STRUCT_T(mouseidx).REGION_ID)

        BW = imbinarize(double(Mice_STRUCT_T(mouseidx).ATLAS2COMET ==regionidx));
        boundaries = bwboundaries(BW);

        if (~isempty(boundaries))
            for k = 1:1

                b = boundaries{k};
                 b(:,2) = b(:,2) + 985 - Bregma_Position(mouseidx,1);
                plot(b(:,2), b(:,1), 'Color' , [0 0 0], 'LineWidth', 1); % Plot boundary with green line
                %     fill(b(:,2)+1, b(:,1)+1, Mice_STRUCT_T(1).COLOR(11,:),'FaceAlpha',0.5,'LineStyle',':');
                hold on

            end
            hold on
        end
        hold on
        regionidx
    end
% 
    hold on
    scatter( Mice_STRUCT_T(mouseidx).COORDINATES(:,1) + 985 - Bregma_Position(mouseidx,1), Mice_STRUCT_T(mouseidx).COORDINATES(:,2),4, Degree_centrality_Cells_color(:,:,learnidx),'filled', 'MarkerFaceAlpha', 0.5 )

hold on
% %pause
end


learnidx

end







% 
% mean_degree_centrality_9=[]
% mean_degree_centrality_10=[]
% mean_degree_centrality_T=[]
% 
% SSp_ul_index = [7   8]
% SSp_ll_index = [21 22]
% SSp_un_index = [23 24]
% SSp_tr_index = [19 20]
% 
% RSPagl_index = [13 14];
% MOs_index    = [11 12];
% MOp_index    = [3   4];
% 
% 
% RSPd_index= [5 6];
% 
% SSP_bfd_index= [17 18];
% VISpm_index= [15 16];
% VISa_index= [25 26];
% 
% VISp_index= [9 10];



%%=========================================================================
%%=========================================================================
%%==================================degrees between cortical regions using
% the neurons activity=======


close all


Cortex_Incommon_New       = Cortex_Incommon([3:4, 1:2, 15:16,  11:12,  13:14, 17:18,9:10, 19:20,    21:22, 25:26, 23:24,5:6, 7:8]);
[~, indexInArray2]        = ismember(Cortex_Incommon_New,Cortex_Incommon_Sorted);


Frame_chosen_index        = Frame_Learn_index;
% Frame_chosen_index        = Frame_CUE_index;



%%-------------START-The total degreee analysis-------------------------

DegreeCentrality_learn_T      = []
DegreeCentrality_cue_T        = []


for mouseidx                  = 1:5

    Trace_Region_T=[]

    for  regionidx            = 1: 26


        Trace_Region10        = [];

        Trace_Region10(:,:)   = Mice_STRUCT_T(mouseidx).NEURON.Trace_unify(Mice_STRUCT_T(mouseidx).Neuron_REGION_ID==CortexID_Mouse_Sorted(mouseidx,indexInArray2(regionidx)) , logical(Frame_Inuse_Mouse{mouseidx}));

        Trace_Region_T        = [Trace_Region_T;Trace_Region10];

    end

    for learnidx              = 1:max(Frame_Learn_index)

        Trace_Region_T_learn  = Trace_Region_T(:,Frame_Learn_index==learnidx);

        CORRE_LEARN_T         = corr(Trace_Region_T_learn');
        CORRE_LEARN_T(isnan(CORRE_LEARN_T))    = 0;
        CORRE_LEARN_T(CORRE_LEARN_T<0)         = 0;

        CORRE_LEARN_T_adjcent = CORRE_LEARN_T - diag(diag(CORRE_LEARN_T));

        DegreeCentrality_learn_T(learnidx,mouseidx) = mean(  sum(CORRE_LEARN_T_adjcent,2)/(size(CORRE_LEARN_T_adjcent,1)-1)  );

        learnidx

    end


    for cueidx                = 1:max(Frame_CUE_index)

        Trace_Region_T_cue    = Trace_Region_T(:,Frame_CUE_index==cueidx);

        CORRE_CUE_T           = corr(Trace_Region_T_cue');
        CORRE_CUE_T(isnan(CORRE_CUE_T))     = 0;
        CORRE_CUE_T(CORRE_CUE_T<0)          = 0;

        CORRE_CUE_T_adjcent   = CORRE_CUE_T - diag(diag(CORRE_CUE_T));

        DegreeCentrality_cue_T(cueidx,mouseidx) = mean(  sum(CORRE_CUE_T_adjcent,2)/(size(CORRE_CUE_T_adjcent,1)-1)  );

        learnidx

    end


    mouseidx



end



%%-------------END-The total degreee analysis------------------------------

%%-------------Cortex regions degree analysis------------------------------

Frame_CUE_index = zeros( 1 , length(Trace_unify_AVG_Norm) );

cuelength = 60

Frame_CUE_index(1,miniscopecuebeginFrame(1):miniscopecuebeginFrame(1)+cuelength) = 1;
Frame_CUE_index(1,miniscopecuebeginFrame(2):miniscopecuebeginFrame(2)+cuelength) = 2;
Frame_CUE_index(1,miniscopecuebeginFrame(3):miniscopecuebeginFrame(3)+cuelength) = 3;
Frame_CUE_index(1,miniscopecuebeginFrame(4):miniscopecuebeginFrame(4)+cuelength) = 4;
Frame_CUE_index(1,miniscopecuebeginFrame(5):miniscopecuebeginFrame(5)+cuelength) = 5;        

Cortex_Incommon_New       = Cortex_Incommon([3:4, 1:2, 15:16,  11:12,   17:18,9:10,13:14, 19:20,    21:22, 25:26, 23:24,5:6, 7:8]);
[~, indexInArray2]        = ismember(Cortex_Incommon_New,Cortex_Incommon_Sorted);

DegreeCentrality_learn_Reg       = []
DegreeCentrality_cue_Reg         = []
Degree_learn_Reg  = []
Degree_cue_Reg    = []


for mouseidx                     = 1:5

    for  regionidxL              = 1: 26

        Trace_Region9            = []
        Trace_Region9(:,:)       = Mice_STRUCT_T(mouseidx).NEURON.Trace_unify(Mice_STRUCT_T(mouseidx).Neuron_REGION_ID==CortexID_Mouse_Sorted(mouseidx,indexInArray2(regionidxL)) , logical(Frame_Inuse_Mouse{mouseidx}));

        for regionidx            = 1: 26

            Trace_Region10       = []

            Trace_Region10(:,:)  = Mice_STRUCT_T(mouseidx).NEURON.Trace_unify(Mice_STRUCT_T(mouseidx).Neuron_REGION_ID==CortexID_Mouse_Sorted(mouseidx,indexInArray2(regionidx))  , logical(Frame_Inuse_Mouse{mouseidx}));

            Trace_Region_T       = [Trace_Region9;Trace_Region10];

            list_in              = 1:size(Trace_Region9,1);
            list_out             = (1+size(Trace_Region9,1)):size(Trace_Region_T,1);

            for learnidx         = 1:max(Frame_Learn_index)

                Trace_Region_T_learn    = Trace_Region_T(:,Frame_Learn_index==learnidx);

                CORRE_learn_Reg    = corr(Trace_Region_T_learn');
                CORRE_learn_Reg(isnan(CORRE_learn_Reg))    = 0;
                CORRE_learn_Reg(CORRE_learn_Reg<0)         = 0;
                CORRE_learn_Reg_adjcent = CORRE_learn_Reg - diag(diag(CORRE_learn_Reg));

                DegreeCentrality_learn_Reg(regionidxL,regionidx,learnidx,mouseidx) = mean(mean(CORRE_learn_Reg_adjcent(list_in , list_out) , 2));
                Degree_learn_Reg(regionidxL,regionidx,learnidx,mouseidx) = sum(mean(CORRE_learn_Reg_adjcent(list_in , list_out) , 2));


            end

            for cueidx         = 1:max(Frame_CUE_index)

                Trace_Region_T_learn    = Trace_Region_T(:,Frame_CUE_index==cueidx);

                CORRE_cue_Reg    = corr(Trace_Region_T_learn');
                CORRE_cue_Reg(isnan(CORRE_cue_Reg))    = 0;
                CORRE_cue_Reg(CORRE_cue_Reg<0)         = 0;
                CORRE_cue_Reg_adjcent = CORRE_cue_Reg - diag(diag(CORRE_cue_Reg));

                DegreeCentrality_cue_Reg(regionidxL,regionidx,cueidx,mouseidx) = mean(mean(CORRE_cue_Reg_adjcent(list_in , list_out) , 2));
                Degree_cue_Reg(regionidxL,regionidx,cueidx,mouseidx) = sum(mean(CORRE_cue_Reg_adjcent(list_in , list_out) , 2));


            end

        end

    end

end


close all

DegreeCentrality_learn_Reg_MiceAVG      = mean( DegreeCentrality_learn_Reg , 4 )  ;
DegreeCentrality_cue_Reg_MiceAVG        = mean( DegreeCentrality_cue_Reg   , 4 )  ;

Degree_learn_Reg_MiceAVG                = mean( Degree_learn_Reg , 4 )  ;
Degree_cue_Reg_MiceAVG                  = mean( Degree_cue_Reg   , 4 )  ;



DegreeCentrality_learn_Reg_MiceAVG_diag = []
DegreeCentrality_cue_Reg_MiceAVG_diag   = []

for learnidx=1:6

    DegreeCentrality_learn_Reg_MiceAVG_diag(:,learnidx) = diag( DegreeCentrality_learn_Reg_MiceAVG(:,:,learnidx)  );
    temp = DegreeCentrality_learn_Reg_MiceAVG(:,:,learnidx);
    Number_temp = size(temp,1);
    [rows, cols] = size( temp );
    [row_indices, col_indices] = ind2sub([rows, cols], 1:numel(temp));

    % 创建新的 binning 索引
    % 例如，原始行的 1,2 映射到新行的 1，3,4 映射到新行的 2
    binned_row_idx = ceil(row_indices / 2);
    binned_col_idx = ceil(col_indices / 2);

    % 将二维索引转换为 accumarray 需要的线性索引
    subs = [binned_row_idx(:), binned_col_idx(:)];

    % 使用 accumarray 对相同 binning 索引的元素求和
    temp_bin = accumarray(subs, temp(:), [Number_temp/2 Number_temp/2], @sum)/4;
     
    DegreeCentrality_learn_Reg_MiceAVG_13Regions(:,:,learnidx) = temp_bin;
    
end
for cueidx=1:5

    DegreeCentrality_cue_Reg_MiceAVG_diag(:,cueidx) = diag( DegreeCentrality_cue_Reg_MiceAVG(:,:,cueidx)  );

     temp = DegreeCentrality_cue_Reg_MiceAVG(:,:,cueidx);
    Number_temp = size(temp,1);
    [rows, cols] = size( temp );
    [row_indices, col_indices] = ind2sub([rows, cols], 1:numel(temp));

    % 创建新的 binning 索引
    % 例如，原始行的 1,2 映射到新行的 1，3,4 映射到新行的 2
    binned_row_idx = ceil(row_indices / 2);
    binned_col_idx = ceil(col_indices / 2);

    % 将二维索引转换为 accumarray 需要的线性索引
    subs = [binned_row_idx(:), binned_col_idx(:)];

    % 使用 accumarray 对相同 binning 索引的元素求和
    temp_bin = accumarray(subs, temp(:), [Number_temp/2 Number_temp/2], @sum)/4;
     
    DegreeCentrality_cue_Reg_MiceAVG_13Regions(:,:,cueidx) = temp_bin;
    
    
end





DegreeCentrality_learn_Reg_MiceAVG_diag_13Reg = ( DegreeCentrality_learn_Reg_MiceAVG_diag(1:2:end,:) + DegreeCentrality_learn_Reg_MiceAVG_diag(2:2:end,:) ) /2;
DegreeCentrality_cue_Reg_MiceAVG_diag_13Reg = ( DegreeCentrality_cue_Reg_MiceAVG_diag(1:2:end,:) + DegreeCentrality_cue_Reg_MiceAVG_diag(2:2:end,:) ) /2;



[DegreeCentrality_learn_Reg_MiceAVG_diag_13Reg_sort index]= sortrows(DegreeCentrality_learn_Reg_MiceAVG_diag_13Reg, 3, 'descend');


Cortex_13regions = Cortex_Incommon_New(1:2:26);
for regionidx=1:13

    temp  = string(  Cortex_13regions(regionidx)  )
    temp2 = extractBetween(temp, 1, strlength(temp) -3   );
    Cortex_13regions_new(regionidx) = string(temp2)
end


Cortex_13regions_new_cell = cellstr(Cortex_13regions_new);

figure(1000)
dataTable = array2table(DegreeCentrality_learn_Reg_MiceAVG_diag_13Reg_sort, ...
    'VariableNames', [{'SPON'}, {'L1'} ,{'L2'} ,{'L3'} ,{'L4'}, {'L5'} ], ...
    'RowNames', [Cortex_13regions_new_cell(index) ]);

h=heatmap(  DegreeCentrality_learn_Reg_MiceAVG_diag_13Reg_sort,'CellLabelFormat', '%.3f')
h.Colormap = turbo; % 更换颜色映射
h.ColorbarVisible = 'off';  % 隐藏 colorbar
set(gcf, 'Position', [100, 100, 600, 500]); % [left, bottom, width, height]
caxis(   [0, 0.18]   ); % 或者使用 clim([0, 1]);
% imagesc(  DegreeCentrality_learn_Reg_MiceAVG_diag_13Reg_sort)
% axis square ;
% xticks(1:6 );
% yticks(1:13);
% xticklabels(  [{'SPON'}, {'L1'} ,{'L2'} ,{'L3'} ,{'L4'}, {'L5'} ]  );
% yticklabels(  [Cortex_13regions_new_cell(index) ]  );
% caxis(   [0, 0.18]   ); % 或者使用 clim([0, 1]);
% colorbar


figure(1001)
dataTable = array2table(DegreeCentrality_learn_Reg_MiceAVG_diag_13Reg_sort, ...
    'VariableNames', [{'SPON'}, {'L1'} ,{'L2'} ,{'L3'} ,{'L4'}, {'L5'} ], ...
    'RowNames', [Cortex_13regions_new_cell(index) ]);

h=heatmap(  DegreeCentrality_learn_Reg_MiceAVG_diag_13Reg,'CellLabelFormat', '%.3f')
h.Colormap = turbo; % 更换颜色映射
h.ColorbarVisible = 'off';  % 隐藏 colorbar
h.FontSize = 14; 
set(gcf, 'Position', [100, 100, 600, 500]); % [left, bottom, width, height]
caxis(   [0, 0.18]   ); % 或者使用 clim([0, 1]);
% imagesc(  DegreeCentrality_learn_Reg_MiceAVG_diag_13Reg)
% axis square ;
% xticks(1:6 );
% yticks(1:13);
% xticklabels(  [{'SPON'}, {'L1'} ,{'L2'} ,{'L3'} ,{'L4'}, {'L5'} ]  );
% yticklabels(  [Cortex_13regions_new_cell ]  );
% caxis(   [0, 0.18]   ); % 或者使用 clim([0, 1]);
% colorbar





figure(1002)
dataTable = array2table(DegreeCentrality_cue_Reg_MiceAVG_diag(index,:), ...
    'VariableNames',  {'C1','C2','C3','C4','C5'}, ...
    'RowNames', [string(Cortex_13regions_new_cell(index))] );

h=heatmap(  DegreeCentrality_cue_Reg_MiceAVG_diag(index,:),'CellLabelFormat', '%.3f')
h.Colormap = turbo; % 更换颜色映射
h.ColorbarVisible = 'off';  % 隐藏 colorbar
h.FontSize = 14; 
set(gcf, 'Position', [100, 100, 500, 500]); % [left, bottom, width, height]
caxis(   [0, 0.18]   ); % 或者使用 clim([0, 1]);
% imagesc(  DegreeCentrality_cue_Reg_MiceAVG_diag(index,:))
% axis square ;
% xticks(1:5 );
% yticks(1:13);
% xticklabels(  [ {'C1'} ,{'C2'} ,{'C3'} ,{'C4'}, {'C5'} ]  );
% yticklabels(  [Cortex_13regions_new_cell(index) ]  );
% caxis(   [0, 0.18]   ); % 或者使用 clim([0, 1]);
% colorbar

figure(1003)
% dataTable = array2table(DegreeCentrality_cue_Reg_MiceAVG_diag, ...
%     'VariableNames',  {'C1','C2','C3','C4','C5'}, ...
%     'RowNames', [string(Cortex_13regions_new_cell(index))] );

h=heatmap(  DegreeCentrality_cue_Reg_MiceAVG_diag_13Reg,'CellLabelFormat', '%.3f')
h.Colormap = turbo; % 更换颜色映射
h.ColorbarVisible = 'off';  % 隐藏 colorbar
h.FontSize = 14; 
set(gcf, 'Position', [100, 100, 500, 500]); % [left, bottom, width, height]
caxis(   [0, 0.2]   ); % 或者使用 clim([0, 1]);
% imagesc(  DegreeCentrality_cue_Reg_MiceAVG_diag(index,:))
% axis square ;
% xticks(1:5 );
% yticks(1:13);
% xticklabels(  [ {'C1'} ,{'C2'} ,{'C3'} ,{'C4'}, {'C5'} ]  );
% yticklabels(  [Cortex_13regions_new_cell(index) ]  );
% caxis(   [0, 0.18]   ); % 或者使用 clim([0, 1]);
% colorbar








for learnidx=1:6


    figure(2000+learnidx)
    h=heatmap(   DegreeCentrality_learn_Reg_MiceAVG_13Regions(:,:,learnidx),'CellLabelColor', 'none')

    disp(['Type of ax: ' class(ax)]);
    h.Colormap = jet ; % 更换颜色映射
    h.ColorbarVisible = 'off';  % 隐藏 colorbar
    h.FontSize = 14;
    h.GridVisible= 'off' ;
    set(gcf, 'Position', [100, 100, 500, 500]); % [left, bottom, width, height]
    caxis(   [0, 0.18]   ); % 或者使用 clim([0, 1]);

    %     ax.GridAlpha = 0.7;       % Slightly transparent
    %     ax.LineWidth = 1;         % Default line width
%       axis off;    % Removes the outer box 
set(ax,'Box', 'off');
end

for learnidx=1:6


    figure(2000+learnidx)
    imagesc(   DegreeCentrality_learn_Reg_MiceAVG_13Regions(:,:,learnidx))
    colormap(jet)
    set(gcf, 'Position', [100, 100, 500, 500]); % [left, bottom, width, height]
    caxis(   [0, 0.18]   ); % 或者使用 clim([0, 1]);
    box off
    axis off
    axis equal

end

 
close all

for learnidx  = 1:6

    temp      = double(DegreeCentrality_learn_Reg_MiceAVG_13Regions(:,:,learnidx));%>mean(DegreeCentrality_learn_Reg_MiceAVG_13Regions(:,:,learnidx),'all'));
    temp_triu = triu(temp);
    temp_1D   = temp_triu(temp_triu ~= 0);

%     figure(2100)
%     histogram(temp_1D)
%     
% 
%     threshold = mean(temp_1D);
%     temp_thr = 1*double(temp>=threshold);
% 
%     figure(2200+learnidx)
%     h=heatmap( temp_thr ,'CellLabelFormat', '%.3f')
%     h.CellLabelColor = 'none'
%     h.Colormap = parula; % 更换颜色映射
%     h.ColorbarVisible = 'off';  % 隐藏 colorbar
%     h.FontSize = 14;
%     set(gcf, 'Position', [100, 100, 500, 500]); % [left, bottom, width, height]
%     caxis(   [0, 1]   ); % 或者使用 clim([0, 1]);


    figure(2300+learnidx)
    h=heatmap( temp ,'CellLabelFormat', '%.3f')
%     h.CellLabelFormat = @(temp) iff(temp >threshold, sprintf('%.3f', x), '');
    h.Colormap = jet; % 更换颜色映射
    h.ColorbarVisible = 'off';  % 隐藏 colorbar
    h.FontSize = 14;
    set(gcf, 'Position', [100, 100, 500, 500]); % [left, bottom, width, height]
    caxis(   [0, 0.18]   ); % 或者使用 clim([0, 1]);
    

end



nodeIds        = [1; 2;3;4;5;6;7;8;9;10;11;12;13]; % 节点 ID

nodeIds_chosen = [2:8]; % 节点 ID
nodeIds_unchosen = nodeIds(~ismember(nodeIds, nodeIds_chosen)); % 节点 ID

nodeLabels     = Cortex_13regions_new_cell(nodeIds); % 节点标签
nodeCategories = {'Group1'; 'Group2'; 'Group2'; 'Group2'; 'Group2'; 'Group2'; 'Group2'; 'Group2';'Group3';'Group3';'Group3';'Group4';'Group5'}; % 节点类别
nodeCategories_chosen = nodeCategories(nodeIds)

color_inuse_GROUP                  = []
color_inuse_GROUP                  = color_inuse([19 3 2 16 24],:);
color_inuse_GROUP(6,:)             = [0.5 0.5 0.5];


temp                               = double(DegreeCentrality_learn_Reg_MiceAVG_13Regions);%>mean(DegreeCentrality_learn_Reg_MiceAVG_13Regions(:,:,learnidx),'all'));
correlation_matrix_T_Mice_Learn_avg_GRAPH = []

temp_rescale= rescale(temp);

close all

for learnidx = 1:6


    temp2 = temp(:,:,learnidx) - diag(diag(temp(:,:,learnidx)));
    temp2(nodeIds_unchosen,:) = 0;


    nodeDegrees = (diag(temp_rescale(:,:,learnidx)));


    temp2(  temp2<1.03*mean(temp,"all")  )         = 0;%1.08*mean(temp,"all")

    correlation_matrix_T_Mice_Learn_avg_GRAPH(:,:,learnidx)  = (temp2);

    [source, target, weight] = find(triu( correlation_matrix_T_Mice_Learn_avg_GRAPH(:,:,learnidx) ));

    edgeList = [source, target, weight];

    G = graph( source, target, weight, Cortex_13regions_new_cell );

    y = [ 6  5    4    3.5  4  2.5  2.5    2     1.2  1  1  2  2]; % 节点的 x 坐标
    x = [ 3  2.8  4    3.5  5  5    3.5    2.8   3.2  2  4  1  5]; % 节点的 y 坐标


    weight2     = rescale(weight);
   
    nodesize    = floor((nodeDegrees)*20+1);
    % 可视化图
   
    handle3=figure(1000+learnidx);
    handle3.Position = [100+learnidx*250 100 500 600];
    hold on;

    nodes = [x;y]'
    % 绘制曲线边
    for i = 1:size(source, 1)
        % 获取边的起点和终点

        if(  ismember( source(i),nodeIds_chosen )  )
           if(ismember( target(i),nodeIds_chosen ))
             startNode = nodes(source(i), :);
             endNode   = nodes(target(i), :);
       
        % 计算控制点（用于贝塞尔曲线）
        controlPoint = (startNode + endNode) / 2 + [0, 0.5];  % 调整控制点位置

        % 生成贝塞尔曲线
        t = linspace(0, 1, 100);
        t = t';

        curve = (1 - t).^2 .* repmat(startNode,length(t),1) + 2 * (1 - t) .* t .* repmat(controlPoint,length(t),1) + t.^2 .* repmat(endNode,length(t),1);
        % 绘制曲线

        plot(curve(:, 1), curve(:, 2), 'k', 'LineWidth', (weight2(i)+0.01)*3  );
          end
         end

    end

    hold on;
    scatter(x(nodeIds_chosen), y(nodeIds_chosen), nodesize(nodeIds_chosen)*15, repmat( color_inuse_GROUP(4,:) , length(nodeIds_chosen),1),'filled');

    text(x(nodeIds_chosen), y(nodeIds_chosen), G.Nodes.Name(nodeIds_chosen), 'VerticalAlignment','bottom', 'HorizontalAlignment','right');
    xlim([0 6])
    ylim([0 7])
    set(gca, 'XColor', 'none');
    set(gca, 'XTickLabel', []);
    set(gca, 'YColor', 'none');
    set(gca, 'YTickLabel', []);

end








for learnidx=1:5


    figure(3500+learnidx)
    h=heatmap(  DegreeCentrality_cue_Reg_MiceAVG_13Regions(:,:,learnidx),'CellLabelFormat', '%.3f')
    h.Colormap = jet; % 更换颜色映射
    h.ColorbarVisible = 'off';  % 隐藏 colorbar
    h.FontSize = 14;
    set(gcf, 'Position', [100, 100, 500, 500]); % [left, bottom, width, height]
    caxis(   [0, 0.35]   ); % 或者使用 clim([0, 1]);


end




close all
for learnidx=1:5

    temp      = double(DegreeCentrality_cue_Reg_MiceAVG_13Regions(:,:,learnidx));%>mean(DegreeCentrality_learn_Reg_MiceAVG_13Regions(:,:,learnidx),'all'));
    temp_triu = triu(temp);
    temp_1D   = temp_triu(temp_triu ~= 0);

    figure(3100)
    histogram(temp_1D)
    

    threshold = mean(temp_1D)*1.02;
    temp_thr = 1*double(temp>=threshold);

    figure(3200+learnidx)
    h=heatmap( temp_thr ,'CellLabelFormat', '%.3f')
    h.CellLabelColor = 'none'
    h.Colormap = parula; % 更换颜色映射
    h.ColorbarVisible = 'off';  % 隐藏 colorbar
    h.FontSize = 14;
    set(gcf, 'Position', [100, 100, 500, 500]); % [left, bottom, width, height]
    caxis(   [0, 1]   ); % 或者使用 clim([0, 1]);


    figure(3300+learnidx)
    h=heatmap( temp ,'CellLabelFormat', '%.3f')
%     h.CellLabelFormat = @(temp) iff(temp >threshold, sprintf('%.3f', x), '');
    h.Colormap = jet; % 更换颜色映射
    h.ColorbarVisible = 'off';  % 隐藏 colorbar
    h.FontSize = 14;
    set(gcf, 'Position', [100, 100, 500, 500]); % [left, bottom, width, height]
    caxis(   [0, 0.35]   ); % 或者使用 clim([0, 1]);
    

end






Frame_CUE_index = zeros( 1 , length(Trace_unify_AVG_Norm) );

cuelength = 60

Frame_CUE_index(1,miniscopecuebeginFrame(1):miniscopecuebeginFrame(1)+cuelength) = 1;
Frame_CUE_index(1,miniscopecuebeginFrame(2):miniscopecuebeginFrame(2)+cuelength) = 2;
Frame_CUE_index(1,miniscopecuebeginFrame(3):miniscopecuebeginFrame(3)+cuelength) = 3;
Frame_CUE_index(1,miniscopecuebeginFrame(4):miniscopecuebeginFrame(4)+cuelength) = 4;
Frame_CUE_index(1,miniscopecuebeginFrame(5):miniscopecuebeginFrame(5)+cuelength) = 5;        

Cortex_Incommon_New       = Cortex_Incommon([3:4, 1:2, 15:16,  11:12,   17:18,9:10,13:14, 19:20,    21:22, 25:26, 23:24,5:6, 7:8]);
[~, indexInArray2]        = ismember(Cortex_Incommon_New,Cortex_Incommon_Sorted);




Frame_chosen_index        = Frame_Learn_index;
% Frame_chosen_index        = Frame_CUE_index

temp_learn_matrix               = []
temp_learn_matrix_T=[]
degree_centrality_out=[]
for mouseidx                    = 1:5

temp                            = []

    for  regionidx              = 1: 26

        temp(regionidx,:)       = mean( Mice_STRUCT_T(mouseidx).NEURON.Trace_unify(Mice_STRUCT_T(mouseidx).Neuron_REGION_ID==CortexID_Mouse_Sorted(mouseidx,indexInArray2(regionidx)) , logical(Frame_Inuse_Mouse{mouseidx}))  , 1 );

    end


    temp_13region =  (  temp(1:2:end,:) + temp(2:2:end,:)  ) /2;


    for learnidx                = 1:max(Frame_chosen_index)
 
        temp_learn              = temp_13region(:,Frame_chosen_index==learnidx);

        temp_learn_matrix       = corr(temp_learn');

        temp_learn_matrix(isnan(temp_learn_matrix))    = 0;
        temp_learn_matrix(temp_learn_matrix<0)         = 0;

        temp_learn_matrix_T(:,:,learnidx, mouseidx)    = temp_learn_matrix;

        temp_learn_matrix_adjcent = temp_learn_matrix - diag(diag(temp_learn_matrix));

        degree_centrality_out(:,learnidx,mouseidx) = (  sum(temp_learn_matrix_adjcent , 2)/( size(temp_learn_matrix_adjcent,1)-1 ) );


    end



end


degree_centrality_out_avgmice = mean(degree_centrality_out,3);

temp_learn_matrix_avg = mean(temp_learn_matrix_T , 4);

close all

for learnidx=1:size(temp_learn_matrix_avg,3)

    figure(4000+learnidx)
    imagesc(temp_learn_matrix_avg(:,:,learnidx))
    xticks(1:13);
    yticks(1:13);
    xticklabels(  [Cortex_13regions_new_cell]  );
    yticklabels(  [Cortex_13regions_new_cell]  );
    set(gca, 'XTick', [], 'YTick', []); % 隐藏刻度
    box off;

    caxis([0, 1]);
    axis off; % 关闭坐标轴和默认 Box
    axis equal;      % 使x和y方向等比例
    axis tight;      % 去除多余空白（可选）
    %     figure(5000+learnidx)
    %     imagesc(  (temp_learn_matrix_avg(1:2:26,1:2:26,learnidx)+temp_learn_matrix_avg(2:2:26,2:2:26,learnidx))/2  )


end


for learnidx=1:size(temp_learn_matrix_avg,3)

    figure(4500+learnidx)
%     h=heatmap( double(temp_learn_matrix_avg(:,:,learnidx)) ,'CellLabelFormat', '%.3f')
%     h.CellLabelColor = 'none'
%     h.Colormap = parula; % 更换颜色映射
%     h.ColorbarVisible = 'off';  % 隐藏 colorbar
%     h.FontSize = 14;
%     set(gcf, 'Position', [100, 100, 500, 500]); % [left, bottom, width, height]
%     caxis(   [0, 1]   ); % 或者使用 clim([0, 1]);
% % axis equal;      % 使x和y方向等比例
% % axis tight;      % 去除多余空白（可选）
% ax = gca;
% % 设置坐标轴位置为正方形
% ax.Position(3) = ax.Position(4); % 使宽度=高度
imagesc(double(temp_learn_matrix_avg(:,:,learnidx)));
colormap(parula);  % 设置颜色映射
        % 添加颜色条

% 启用并自定义网格

axis equal;      % 使x和y方向等比例
axis tight;      % 去除多余空白（可选）
caxis(   [0, 1]   );
set(gca, 'XTick', [], 'YTick', []); % 隐藏刻度
% [x, y] = meshgrid(1:size(temp_learn_matrix_avg, 2)+1, 1:size(temp_learn_matrix_avg, 1)+1);
% 
% % 绘制自定义颜色的网格线
% hold on;
% for i = 1:size(temp_learn_matrix_avg, 1)
%     plot([0.5, size(temp_learn_matrix_avg, 2)+0.5], [i+0.5, i+0.5], 'w-', 'LineWidth', .1); % 红色横线
% end
% for j = 1:size(temp_learn_matrix_avg, 2)
%     plot([j+0.5, j+0.5], [0.5, size(temp_learn_matrix_avg, 1)+0.5], 'w-', 'LineWidth', .1); % 蓝色竖线
% end
% hold off;
end




nodeIds        = [1; 2; 3; 4;5;6;7;8;9;10;11;12;13]; % 节点 ID
nodeLabels     = Cortex_13regions_new_cell; % 节点标签
nodeCategories = {'Group1'; 'Group2'; 'Group2'; 'Group2'; 'Group2'; 'Group2'; 'Group2'; 'Group2';'Group3';'Group3';'Group3';'Group4';'Group5'}; % 节点类别

color_inuse_GROUP                  = []
color_inuse_GROUP                  = color_inuse([19 3 2 16 24],:);
color_inuse_GROUP(6,:)             = [0.5 0.5 0.5];


temp                               = double(DegreeCentrality_learn_Reg_MiceAVG_13Regions);%>mean(DegreeCentrality_learn_Reg_MiceAVG_13Regions(:,:,learnidx),'all'));


close all
for learnidx = 1:6

    temp2 = temp(:,:,learnidx) - diag(diag(temp(:,:,learnidx)));
    
    nodeDegrees = (diag(temp_rescale(:,:,learnidx)));


%     temp2(  temp2<prctile(reshape(temp(:,:,learnidx),size(temp(:,:,learnidx),1)*size(temp(:,:,learnidx),2),1), 95)  )         = 0;%1.08*mean(temp,"all")%0.046
    temp2(  temp2<max(max(temp2))*0.9 )         = 0;%1.08*mean(temp,"all")%0.046
    correlation_matrix_T_Mice_Learn_avg_GRAPH(:,:,learnidx)  = (temp2);

    [source, target, weight] = find(triu( correlation_matrix_T_Mice_Learn_avg_GRAPH(:,:,learnidx) ));

    edgeList = [source, target, weight];

    G = graph( source, target, weight, Cortex_13regions_new_cell );


    y = [ 5.5  4.5  4.0  3.5  4.0  3.0  2.5  1.5   1.0  1.0  0.5  2.0  3.0]; % 节点的 x 坐标
    x = [ 2.8  3.5  4.5  4.0  5.5  5.5  4.0  3.5   2.8  1.5  2.5  1.0  1.0]; % 节点的 y 坐标


    weight2 =rescale(weight);
   
    nodesize    = floor((nodeDegrees)*10+1);
    % 可视化图
   
    handle3=figure(900+learnidx);
    handle3.Position = [100+learnidx*250 100 500 600];
    hold on;

    nodes = [x;y]'
    % 绘制曲线边
    for i = 1:size(source, 1)
        % 获取边的起点和终点
        startNode = nodes(source(i), :);
        endNode   = nodes(target(i), :);

        % 计算控制点（用于贝塞尔曲线）
        controlPoint = (startNode + endNode) / 2 + [0, 0.5];  % 调整控制点位置

        % 生成贝塞尔曲线
        t = linspace(0, 1, 20);
        t = t';

        curve = (1 - t).^2 .* repmat(startNode,length(t),1) + 2 * (1 - t) .* t .* repmat(controlPoint,length(t),1) + t.^2 .* repmat(endNode,length(t),1);
        % 绘制曲线

%         plot(curve(:, 1), curve(:, 2), 'k.', 'LineWidth',  0.1 );%(weight2(i)+0.01)*1
 plot(curve(:, 1), curve(:, 2), 'k-', 'LineWidth',  (weight2(i)+0.01)*2 );%(weight2(i)+0.01)*1
    end

    hold on;
 

%  rectangle('Position', [2.8-0.75, 5.5-0.25, 1, 0.5], ...
%            'FaceColor', 'none', ...
%            'EdgeColor', color_inuse_GROUP(1,:), ...
%            'LineWidth', 2, ...
%            'Curvature', [0.2, 0.2]); % 圆角矩形
% 
%  hold on;
%  
% 
%     rectangle('Position', [3.2-0.25, 2.5-0.25, 3.0, 2.5], ...
%            'FaceColor', 'none', ...
%            'EdgeColor', color_inuse_GROUP(4,:), ...
%            'LineWidth', 2, ...
%            'Curvature', [0.2, 0.2]); % 圆角矩形
%      hold on;
%  
% 
%     rectangle('Position', [1.0-0.75, 2.0-0.25, 1, 0.5], ...
%            'FaceColor', 'none', ...
%            'EdgeColor', color_inuse_GROUP(2,:), ...
%            'LineWidth', 2, ...
%            'Curvature', [0.2, 0.2]); % 圆角矩形
%      hold on;
%  
% 
%     rectangle('Position', [1.0-0.75, 3.0-0.25, 1, 0.5], ...
%            'FaceColor', 'none', ...
%            'EdgeColor', color_inuse_GROUP(5,:), ...
%            'LineWidth', 2, ...
%            'Curvature', [0.2, 0.2]); % 圆角矩形
%      hold on;
%  
% 
%     rectangle('Position', [1.0-0.25, 0.5-0.25, 2.5, 1.2], ...
%            'FaceColor', 'none', ...
%            'EdgeColor', color_inuse_GROUP(3,:), ...
%            'LineWidth', 2, ...
%            'Curvature', [0.2, 0.2]); % 圆角矩形
%      hold on;
%  



    
    scatter(x, y, nodesize*15, [color_inuse_GROUP(1,:);
        color_inuse_GROUP(4,:);
        color_inuse_GROUP(4,:);
        color_inuse_GROUP(4,:);
        color_inuse_GROUP(4,:);
        color_inuse_GROUP(4,:);
        color_inuse_GROUP(4,:);
        color_inuse_GROUP(4,:);
        color_inuse_GROUP(3,:);
        color_inuse_GROUP(3,:);
        color_inuse_GROUP(3,:);
        color_inuse_GROUP(2,:);
        color_inuse_GROUP(5,:)
        ] , 'filled');

     %text(x, y, G.Nodes.Name, 'VerticalAlignment','bottom', 'HorizontalAlignment','right');
    xlim([0 6])
    ylim([0 7])
    set(gca, 'XColor', 'none');
    set(gca, 'XTickLabel', []);
    set(gca, 'YColor', 'none');
    set(gca, 'YTickLabel', []);

   



end


%%=========================================================================
%%=====================================
%%====================


correlation_matrix_T_Mice_Learn_avg_GRAPH = []

temp_rescale         = rescale(temp);

weight4              = []
weight3              = []
weight2              = []
weight1              = []

Degree_REGION5_learn = []


for learnidx                  = 1:6



    temp2                     = temp(:,:,learnidx) - diag(diag(temp(:,:,learnidx)));
    nodeIds_chosen            = [2:8];
    nodeIds_unchosen          = nodeIds(~ismember(nodeIds, nodeIds_chosen)); % 节点 ID
    temp2(nodeIds_unchosen,:) = 0;
%     temp2(:,nodeIds_unchosen) = 0;

    
    %     temp2(  temp2<0.5*mean(temp,"all")  )         = 0;%1.08*mean(temp,"all")%0.046

    [source, target, weight]  = find(triu( temp2 ));
    weight1(:,learnidx)       = (weight);


    temp2                     = temp(:,:,learnidx) - diag(diag(temp(:,:,learnidx)));
    nodeIds_chosen            = [8:11];
    nodeIds_unchosen          = nodeIds(~ismember(nodeIds, nodeIds_chosen)); % 节点 ID
    temp2(nodeIds_unchosen,:) = 0;
%     temp2(:,nodeIds_unchosen) = 0;

    %     temp2(  temp2<0.5*mean(temp,"all")  )         = 0;

    [source, target, weight]  = find(triu( temp2 ));
    weight2(:,learnidx)       = (weight);


%     temp2                     = temp(:,:,learnidx) - diag(diag(temp(:,:,learnidx)));
%     nodeIds_chosen            = [7,8];
%     nodeIds_unchosen          = nodeIds(~ismember(nodeIds, nodeIds_chosen)); % 节点 ID
%     temp2(nodeIds_unchosen,:) = 0;
%     temp2(:,nodeIds_unchosen) = 0;
% 
%     %     temp2(  temp2<0.5*mean(temp,"all")  )         = 0;
% 
%     [source, target, weight]  = find(triu( temp2 ));
%     weight3(:,learnidx)       = (weight);
% 
%     temp2                     = temp(:,:,learnidx) - diag(diag(temp(:,:,learnidx)));
%     nodeIds_chosen            = [7,9];
%     nodeIds_unchosen          = nodeIds(~ismember(nodeIds, nodeIds_chosen)); % 节点 ID
%     temp2(nodeIds_unchosen,:) = 0;
%     temp2(:,nodeIds_unchosen) = 0;
% 
%     %     temp2(  temp2<0.5*mean(temp,"all")  )         = 0;
% 
%     [source, target, weight]  = find(triu( temp2 ));
%     weight4(:,learnidx)       = (weight);



    Degree_MOs_RSPd           = mean( temp(1,    13,learnidx)  , 'all');
    Degree_MOs_RSPagl         = mean( temp(1,    12,learnidx)  , 'all');
    Degree_MOs_SSP            = mean( temp(1,  2: 7,learnidx)  , 'all');

   
    Degree_MOs_VISa           = mean( temp(1,     8,learnidx)  , 'all');
    Degree_MOs_VISam          = mean( temp(1,     9,learnidx)  , 'all');    
    Degree_MOs_VIS            = mean( temp(1,  8:11,learnidx)  , 'all');


    Degree_SSP_RSPd           = mean( temp(2:7,   13,learnidx) , 'all');
    Degree_SSP_RSPagl         = mean( temp(2:7,   12,learnidx) , 'all');


    Degree_SSP_VISa           = mean( temp(2:7,    8,learnidx) , 'all');
    Degree_SSP_VISam          = mean( temp(2:7,    9,learnidx) , 'all');
    Degree_SSP_VIS            = mean( temp(2:7, 8:11,learnidx) , 'all');


    Degree_RSPd_RSPagl        = mean( temp(13,    12,learnidx) , 'all');

    Degree_RSPd_VIS           = mean( temp(13,  8:11,learnidx) , 'all');
  


    Degree_RSPagl_VISa        = mean( temp(12,     8,learnidx) , 'all');
    Degree_RSPagl_VISam       = mean( temp(12,     9,learnidx) , 'all');
    Degree_RSPagl_VIS         = mean( temp(12,  8:11,learnidx) , 'all');

    

    Degree_REGION5     = ( [ Degree_MOs_SSP     Degree_MOs_RSPd     Degree_MOs_RSPagl     Degree_MOs_VISa Degree_MOs_VISam  Degree_MOs_VIS ...
                             Degree_SSP_VISa    Degree_SSP_VISam    Degree_SSP_VIS        Degree_SSP_RSPd Degree_SSP_RSPagl  ...
                             Degree_RSPd_RSPagl Degree_RSPd_VISa    Degree_RSPd_VISam     Degree_RSPd_VIS ...
                             Degree_RSPagl_VISa Degree_RSPagl_VISam Degree_RSPagl_VIS ] );


    Degree_REGION5_learn(:,learnidx) = Degree_REGION5;



end



Degree_REGION5_Total_norm        = []
Degree_REGION5_Total             = [weight1;weight2;weight3;weight4;Degree_REGION5_learn]
Degree_REGION5_Total2            = Degree_REGION5_Total.^7;
TEMP                             = Degree_REGION5_Total2(:,2:6);
Degree_REGION5_Total_norm(:,1)   = 4*rescale(Degree_REGION5_Total2(:,1));
Degree_REGION5_Total_norm(:,2:6) = 4*(TEMP - min(TEMP(:))) / (max(TEMP(:)) - min(TEMP(:)));

% Degree_REGION5_Total_norm(Degree_REGION5_Total_norm<=0.3) = nan;




close all




COLOR_spon  = [0.2 0.2 0.2]/1
COLOR_learn = [0 0 1]/1

for learnidx   = 1:6

    nodeDegrees = (diag(temp_rescale(:,:,learnidx)));

    x = [ 2.3  3.5  4.5  4.0  5.5  5.5  4.2  3.5   2.5  2.5  3.5  1.0  1.0]; % 节点的 y 坐标
    y = [ 5.5  4.5  4.0  3.5  4.0  3.0  2.8  1.5   1.5  0.5  0.5  2.0  3.0]; % 节点的 x 坐标

    x(2:7) = x(2:7) +0.3;

   if learnidx>=4

        x(8) = 3.8
        y(8) = 2.5
   end


   if learnidx==1

       base_linewidth = 0.1;
       linestyle = '--'

   else
       base_linewidth = 0.01;
       linestyle = '-'
   end


    temp2                     = temp(:,:,learnidx) - diag(diag(temp(:,:,learnidx)));
    nodeIds_chosen            = [2:8]; % 节点 ID
    nodeIds_unchosen          = nodeIds(~ismember(nodeIds, nodeIds_chosen)); % 节点 ID
    temp2(nodeIds_unchosen,:) = 0;
%     temp2(:,nodeIds_unchosen) = 0;

    correlation_matrix_T_Mice_Learn_avg_GRAPH(:,:,learnidx)  = (temp2);

    [source, target, weight] = find(triu( temp2 ));

    edgeList = [source, target, weight];

    G = graph( source, target, weight, Cortex_13regions_new_cell );



    nodesize    = floor((nodeDegrees)*20+1);


    handle3=figure(1000+learnidx);
    handle3.Position = [100+learnidx*250 100 550 600];
    hold on;

    nodes = [x;y]'
    % 绘制曲线边
    for i = 1:size(source, 1)
        % 获取边的起点和终点

        if(  ismember( source(i),nodeIds_chosen )  )
            if(ismember( target(i),nodeIds_chosen ))
                startNode = nodes(source(i), :);
                endNode   = nodes(target(i), :);
            end
        end
        % 计算控制点（用于贝塞尔曲线）
        controlPoint = (startNode + endNode) / 2 + [0, 0.5];  % 调整控制点位置

        % 生成贝塞尔曲线
        t = linspace(0, 1, 40);
        t = t';

        curve = (1 - t).^2 .* repmat(startNode,length(t),1) + 2 * (1 - t) .* t .* repmat(controlPoint,length(t),1) + t.^2 .* repmat(endNode,length(t),1);
        % 绘制曲线

        %         plot(curve(:, 1), curve(:, 2), '--','Color',COLOR_spon, 'LineWidth', (Degree_REGION5_Total_norm(i,learnidx)+0.01)*3  );

        plot(curve(:, 1), curve(:, 2), linestyle,'Color',COLOR_spon, 'LineWidth', (Degree_REGION5_Total_norm(i,learnidx)+base_linewidth)  );

    end


    hold on


    temp2 = temp(:,:,learnidx) - diag(diag(temp(:,:,learnidx)));
    nodeIds_chosen = [8:11]; % 节点 ID
    nodeIds_unchosen = nodeIds(~ismember(nodeIds, nodeIds_chosen)); % 节点 ID
    temp2(nodeIds_unchosen,:) = 0;
%     temp2(:,nodeIds_unchosen) = 0;
 
%     temp2(  temp2<1.03*mean(temp,"all")  )         = 0;%1.08*mean(temp,"all")
    %     temp2(  temp2<1*mean(temp,"all")  )         = 0;
 
    correlation_matrix_T_Mice_Learn_avg_GRAPH(:,:,learnidx)  = (temp2);

    [source, target, weight] = find(triu( temp2 ));

    edgeList = [source, target, weight];

    G = graph( source, target, weight, Cortex_13regions_new_cell );



    nodesize    = floor((nodeDegrees)*20+1);



    % 绘制曲线边
    for i = 1:size(source, 1)
        % 获取边的起点和终点

        if(  ismember( source(i),nodeIds_chosen )  )
            if(ismember( target(i),nodeIds_chosen ))
                startNode = nodes(source(i), :);
                endNode   = nodes(target(i), :);
            end
        end
        % 计算控制点（用于贝塞尔曲线）
        controlPoint = (startNode + endNode) / 2 + [0, 0.5];  % 调整控制点位置

        % 生成贝塞尔曲线
        t = linspace(0, 1, 40);
        t = t';

        curve = (1 - t).^2 .* repmat(startNode,length(t),1) + 2 * (1 - t) .* t .* repmat(controlPoint,length(t),1) + t.^2 .* repmat(endNode,length(t),1);
        % 绘制曲线

        %         plot(curve(:, 1), curve(:, 2), '--','Color',COLOR_spon, 'LineWidth', (Degree_REGION5_Total_norm(i+length(weight1),learnidx)+0.01)   );
        %
        plot(curve(:, 1), curve(:, 2),linestyle, 'Color',COLOR_spon, 'LineWidth', (Degree_REGION5_Total_norm(i+length(weight1),learnidx)+base_linewidth)   );

    end

    hold on

% 
%     temp2 = temp(:,:,learnidx) - diag(diag(temp(:,:,learnidx)));
%     nodeIds_chosen = [7,8]; % 节点 ID
%     nodeIds_unchosen = nodeIds(~ismember(nodeIds, nodeIds_chosen)); % 节点 ID
%     temp2(nodeIds_unchosen,:) = 0;
%     temp2(:,nodeIds_unchosen) = 0;
%     %     temp2(  temp2<1.03*mean(temp,"all")  )         = 0;%1.08*mean(temp,"all")
%     %     temp2(  temp2<1*mean(temp,"all")  )         = 0;
%     correlation_matrix_T_Mice_Learn_avg_GRAPH(:,:,learnidx)  = (temp2);
% 
%     [source, target, weight] = find(triu( temp2 ));
% 
%     edgeList = [source, target, weight];
% 
%     G = graph( source, target, weight, Cortex_13regions_new_cell );
% 
% 
% 
%     nodesize    = floor((nodeDegrees)*20+1);
% 
% 
% 
%     % 绘制曲线边
%     for i = 1:size(source, 1)
%         % 获取边的起点和终点
% 
%         if(  ismember( source(i),nodeIds_chosen )  )
%             if(ismember( target(i),nodeIds_chosen ))
%                 startNode = nodes(source(i), :);
%                 endNode   = nodes(target(i), :);
%             end
%         end
%         % 计算控制点（用于贝塞尔曲线）
%         controlPoint = (startNode + endNode) / 2 + [0, 0.5];  % 调整控制点位置
% 
%         % 生成贝塞尔曲线
%         t = linspace(0, 1, 40);
%         t = t';
% 
%         curve = (1 - t).^2 .* repmat(startNode,length(t),1) + 2 * (1 - t) .* t .* repmat(controlPoint,length(t),1) + t.^2 .* repmat(endNode,length(t),1);
%         % 绘制曲线
% 
%         %         plot(curve(:, 1), curve(:, 2), '--','Color',COLOR_spon, 'LineWidth', (Degree_REGION5_Total_norm(i+length(weight1)+length(weight2),learnidx)+0.01)   );
% 
%         plot(curve(:, 1), curve(:, 2), 'Color',COLOR_learn, 'LineWidth', (Degree_REGION5_Total_norm(i+length(weight1)+length(weight2),learnidx)+0.01)   );
% 
%     end
% 
% 
% 
%     hold on
% 
%     temp2 = temp(:,:,learnidx) - diag(diag(temp(:,:,learnidx)));
%     nodeIds_chosen = [7,9]; % 节点 ID
%     nodeIds_unchosen = nodeIds(~ismember(nodeIds, nodeIds_chosen)); % 节点 ID
%     temp2(nodeIds_unchosen,:) = 0;
%     temp2(:,nodeIds_unchosen) = 0;
%     %     temp2(  temp2<1.03*mean(temp,"all")  )         = 0;%1.08*mean(temp,"all")
%     %     temp2(  temp2<1*mean(temp,"all")  )         = 0;
%     correlation_matrix_T_Mice_Learn_avg_GRAPH(:,:,learnidx)  = (temp2);
% 
%     [source, target, weight] = find(triu( temp2 ));
% 
%     edgeList = [source, target, weight];
% 
%     G = graph( source, target, weight, Cortex_13regions_new_cell );
% 
% 
% 
%     nodesize    = floor((nodeDegrees)*20+1);
% 
% 
% 
%     % 绘制曲线边
%     for i = 1:size(source, 1)
%         % 获取边的起点和终点
% 
%         if(  ismember( source(i),nodeIds_chosen )  )
%             if(ismember( target(i),nodeIds_chosen ))
%                 startNode = nodes(source(i), :);
%                 endNode   = nodes(target(i), :);
%             end
%         end
%         % 计算控制点（用于贝塞尔曲线）
%         controlPoint = (startNode + endNode) / 2 + [0, 0.5];  % 调整控制点位置
% 
%         % 生成贝塞尔曲线
%         t = linspace(0, 1, 40);
%         t = t';
% 
%         curve = (1 - t).^2 .* repmat(startNode,length(t),1) + 2 * (1 - t) .* t .* repmat(controlPoint,length(t),1) + t.^2 .* repmat(endNode,length(t),1);
%         % 绘制曲线
% 
%         %         plot(curve(:, 1), curve(:, 2), '--','Color',COLOR_spon, 'LineWidth', (Degree_REGION5_Total_norm(i+length(weight1)+length(weight2),learnidx)+0.01)   );
% 
%         plot(curve(:, 1), curve(:, 2), 'Color',COLOR_learn, 'LineWidth', (Degree_REGION5_Total_norm(i+length(weight1)+length(weight2)+length(weight3),learnidx)+0.01)   );
% 
%     end



    hold on




    %     Degree_REGION5     = ( [ Degree_MOs_SSP     Degree_MOs_RSPd     Degree_MOs_RSPagl     Degree_MOs_VISa Degree_MOs_VISam  Degree_MOs_VIS ...
    %                              Degree_SSP_VISa    Degree_SSP_VISam    Degree_SSP_VIS        Degree_SSP_RSPd Degree_SSP_RSPagl  ...
    %                              Degree_RSPd_RSPagl Degree_RSPd_VISa    Degree_RSPd_VISam     Degree_RSPd_VIS ...
    %                              Degree_RSPagl_VISa Degree_RSPagl_VISam Degree_RSPagl_VIS ] );

    LENGTH_nodes = size(weight1,1) + size(weight2,1) + size(weight3,1)+ size(weight4,1);


    plot( [2.55 4.45], [5.5  4.75], linestyle,  'Color',  COLOR_spon , 'LineWidth', (Degree_REGION5_Total_norm(LENGTH_nodes+1,learnidx)+ base_linewidth)   );
    hold on
    plot( [2.05 1.25], [5.25 3],  linestyle,    'Color',  COLOR_spon , 'LineWidth', (Degree_REGION5_Total_norm(LENGTH_nodes+2,learnidx)+ base_linewidth)   );
    hold on
    plot( [2.05 1.25], [5.25 2],   linestyle,   'Color',  COLOR_spon , 'LineWidth', (Degree_REGION5_Total_norm(LENGTH_nodes+3,learnidx)+ base_linewidth)   );
    hold on

    plot( [2.05 x(8)], [5.25 y(8)],  linestyle,   'Color', COLOR_learn , 'LineWidth', (Degree_REGION5_Total_norm(LENGTH_nodes+4,learnidx)+ base_linewidth)   );
    hold on
    plot( [2.05 2.5], [5.25 1.5],  linestyle,   'Color', COLOR_learn , 'LineWidth', (Degree_REGION5_Total_norm(LENGTH_nodes+5,learnidx)+ base_linewidth)   );
    hold on
    plot( [2.05 2.25], [5.25 1.05], linestyle,   'Color', COLOR_spon ,  'LineWidth', (Degree_REGION5_Total_norm(LENGTH_nodes+6,learnidx)+ base_linewidth)   );
    hold on



    plot( [4.75 x(8)], [2.25 y(8)], linestyle,  'Color', COLOR_learn , 'LineWidth', (Degree_REGION5_Total_norm(LENGTH_nodes+8,learnidx)+ base_linewidth)    );
    hold on
    plot( [3.45 2.5], [3.5 1.5], linestyle,  'Color', COLOR_learn , 'LineWidth', (Degree_REGION5_Total_norm(LENGTH_nodes+7,learnidx)+ base_linewidth)    );
    hold on
    plot( [4.75 3.75], [2.25 1.05], linestyle,  'Color', COLOR_spon , 'LineWidth', (Degree_REGION5_Total_norm(LENGTH_nodes+9,learnidx)+ base_linewidth)    );
    hold on

    plot( [3.45 1.25], [3.5 3],   linestyle,    'Color', COLOR_spon, 'LineWidth', (Degree_REGION5_Total_norm(LENGTH_nodes+10,learnidx)+ base_linewidth)     );
    hold on
    plot( [3.45 1.25], [3.5 2],   linestyle,     'Color',COLOR_spon, 'LineWidth', (Degree_REGION5_Total_norm(LENGTH_nodes+11,learnidx)+ base_linewidth)     );
    hold on



    plot( [0.75 0.75], [2.75 2.25],  linestyle, 'Color', COLOR_spon, 'LineWidth', (Degree_REGION5_Total_norm(LENGTH_nodes+12,learnidx)+ base_linewidth)     );
    hold on

    plot( [1.25 2.25], [3    1.05],linestyle,   'Color', COLOR_spon, 'LineWidth', (Degree_REGION5_Total_norm(LENGTH_nodes+15,learnidx)+ base_linewidth)     );
    hold on



    plot( [1.25 x(8)], [2 y(8)],   linestyle,   'Color', COLOR_learn, 'LineWidth', (Degree_REGION5_Total_norm(LENGTH_nodes+16,learnidx)+ base_linewidth)    );
    plot( [1.25 2.5], [2 1.5],    linestyle,  'Color', COLOR_learn, 'LineWidth', (Degree_REGION5_Total_norm(LENGTH_nodes+17,learnidx)+ base_linewidth)    );

    plot( [1.25 2.25], [2 1.05],  linestyle,   'Color', COLOR_spon, 'LineWidth', (Degree_REGION5_Total_norm(LENGTH_nodes+18,learnidx)+ base_linewidth)    );


    hold on




    rectangle('Position', [2.3-0.75, 5.5-0.25, 1, 0.5], ...
        'FaceColor', 'none', ...
        'EdgeColor', color_inuse_GROUP(1,:), ...
        'LineWidth', 3, ...
        'Curvature', [0.2, 0.2], ...
        'LineStyle', '-'); % 圆角矩形

    hold on;


    rectangle('Position', [3.7-0.25, 2.5-0.25, 2.6, 2.5], ...
        'FaceColor', 'none', ...
        'EdgeColor', color_inuse_GROUP(4,:), ...
        'LineWidth', 3, ...
        'Curvature', [0.2, 0.2], ...
        'LineStyle', '-'); % 圆角矩形
    hold on;


    rectangle('Position', [1.0-0.75, 2.0-0.25, 1, 0.5], ...
        'FaceColor', 'none', ...
        'EdgeColor', color_inuse_GROUP(2,:), ...
        'LineWidth', 3, ...
        'Curvature', [0.2, 0.2], ...
        'LineStyle', '-'); % 圆角矩形
    hold on;


    rectangle('Position', [1.0-0.75, 3.0-0.25, 1, 0.5], ...
        'FaceColor', 'none', ...
        'EdgeColor', color_inuse_GROUP(5,:), ...
        'LineWidth', 3, ...
        'Curvature', [0.2, 0.2], ...
        'LineStyle', '-'); % 圆角矩形
    hold on;


    rectangle('Position', [2.5-0.25, 0.5-0.25, 1.5, 1.6], ...
        'FaceColor', 'none', ...
        'EdgeColor', color_inuse_GROUP(3,:), ...
        'LineWidth', 3, ...
        'Curvature', [0.2, 0.2], ...
        'LineStyle', '-'); % 圆角矩形
    hold on;



    scatter(x, y, nodesize*15, [color_inuse_GROUP(1,:);
        color_inuse_GROUP(4,:);
        color_inuse_GROUP(4,:);
        color_inuse_GROUP(4,:);
        color_inuse_GROUP(4,:);
        color_inuse_GROUP(4,:);
        color_inuse_GROUP(4,:);
        color_inuse_GROUP(4,:);
        color_inuse_GROUP(3,:);
        color_inuse_GROUP(3,:);
        color_inuse_GROUP(3,:);
        color_inuse_GROUP(2,:);
        color_inuse_GROUP(5,:)
        ] , 'filled');

    %text(x, y, G.Nodes.Name, 'VerticalAlignment','bottom', 'HorizontalAlignment','right');
    xlim([0 6])
    ylim([0 7])
    set(gca, 'XColor', 'none');
    set(gca, 'XTickLabel', []);
    set(gca, 'YColor', 'none');
    set(gca, 'YTickLabel', []);

end


%%====================   
%%=====================================
%%=========================================================================

correlation_matrix_T_Mice_Learn_avg_GRAPH = []
Degree_REGION5_learn=[]
temp_rescale= rescale(temp);

weight2 = []
weight1=[]


for learnidx = 1:6


    temp2 = temp(:,:,learnidx) - diag(diag(temp(:,:,learnidx)));
    nodeIds_chosen   = [2:7];
    nodeIds_unchosen = nodeIds(~ismember(nodeIds, nodeIds_chosen)); % 节点 ID
    temp2(nodeIds_unchosen,:) = 0;
    temp2(:,nodeIds_unchosen) = 0;
    %     temp2(  temp2<0.5*mean(temp,"all")  )         = 0;%1.08*mean(temp,"all")%0.046
    [source, target, weight] = find(triu( temp2 ));
    weight1(:,learnidx)     = (weight);


    temp2 = temp(:,:,learnidx) - diag(diag(temp(:,:,learnidx)));
    nodeIds_chosen = [8:11];
    nodeIds_unchosen = nodeIds(~ismember(nodeIds, nodeIds_chosen)); % 节点 ID
    temp2(nodeIds_unchosen,:) = 0;
    temp2(:,nodeIds_unchosen) = 0;
    %     temp2(  temp2<0.5*mean(temp,"all")  )         = 0;
    [source, target, weight] = find(triu( temp2 ));
    weight2(:,learnidx)     = (weight);


    Degree_MOs_RSPd                    = mean( temp(1,13,learnidx) , 'all');
    Degree_MOs_RSPagl                  = mean( temp(1,12,learnidx) , 'all');
    Degree_MOs_SSP                     = mean( temp(1,2:7,learnidx) , 'all');
    Degree_MOs_VIS                     = mean( temp(1,8:11,learnidx) , 'all');

    Degree_SSP_RSPd                    = mean( temp(2:7,13,learnidx) , 'all');
    Degree_SSP_RSPagl                  = mean( temp(2:7,12,learnidx) , 'all');
    Degree_SSP_VIS                     = mean( temp(2:7,8:11,learnidx) , 'all');

    Degree_RSPd_RSPagl                  = mean( temp(13,12,learnidx) , 'all');
    Degree_RSPd_VIS                     = mean( temp(13,8:11,learnidx) , 'all');

    Degree_RSPagl_VIS                   = mean( temp(12,8:11,learnidx) , 'all');

    Degree_REGION5     = ([Degree_MOs_SSP Degree_MOs_RSPd Degree_MOs_RSPagl  Degree_MOs_VIS ...
        Degree_SSP_VIS Degree_SSP_RSPd Degree_SSP_RSPagl  ...
        Degree_RSPd_RSPagl Degree_RSPd_VIS ...
        Degree_RSPagl_VIS]);

    Degree_REGION5_learn(:,learnidx) = Degree_REGION5;


end


Degree_REGION5_Total_norm = []
Degree_REGION5_Total = [weight1;weight2;Degree_REGION5_learn]
Degree_REGION5_Total2 = Degree_REGION5_Total.^5;
TEMP = Degree_REGION5_Total2(:,2:6);
Degree_REGION5_Total_norm(:,1) = 2*rescale(Degree_REGION5_Total2(:,1));
Degree_REGION5_Total_norm(:,2:6) = 2*(TEMP - min(TEMP(:))) / (max(TEMP(:)) - min(TEMP(:)));

% Degree_REGION5_Total_norm(Degree_REGION5_Total_norm<=0.3) = nan;



close all


COLOR_spon  = [0 0 0]/255
COLOR_learn = [0 0 0]/255

 learnidx = 1

    nodeDegrees = (diag(temp_rescale(:,:,learnidx)));

    y = [ 5.5  4.5  4.0  3.5  4.0  3.0  3.0  1.0   1.0  1.0  0.5  2.0  3.0]; % 节点的 x 坐标
    x = [ 2.8  3.5  4.5  4.0  5.5  5.5  4.0  3.5   2.8  1.5  2.5  1.0  1.0]; % 节点的 y 坐标


 
    temp2 = temp(:,:,learnidx) - diag(diag(temp(:,:,learnidx)));
    nodeIds_chosen = [2:8]; % 节点 ID
    nodeIds_unchosen = nodeIds(~ismember(nodeIds, nodeIds_chosen)); % 节点 ID
    temp2(nodeIds_unchosen,:) = 0;
     temp2(:,nodeIds_unchosen) = 0;
%     temp2(  temp2<1.03*mean(temp,"all")  )         = 0;%1.08*mean(temp,"all")
%     temp2(  temp2<1*mean(temp,"all")  )         = 0;
    correlation_matrix_T_Mice_Learn_avg_GRAPH(:,:,learnidx)  = (temp2);

    [source, target, weight] = find(triu( temp2 ));

    edgeList = [source, target, weight];

    G = graph( source, target, weight, Cortex_13regions_new_cell );
   
 
   
    nodesize    = floor((nodeDegrees)*20+1);
   
   
    handle3=figure(1000+learnidx);
    handle3.Position = [100+learnidx*250 100 500 600];
    hold on;

    nodes = [x;y]'
    % 绘制曲线边
    for i = 1:size(source, 1)
        % 获取边的起点和终点

        if(  ismember( source(i),nodeIds_chosen )  )
           if(ismember( target(i),nodeIds_chosen ))
             startNode = nodes(source(i), :);
             endNode   = nodes(target(i), :);
           end
         end
        % 计算控制点（用于贝塞尔曲线）
        controlPoint = (startNode + endNode) / 2 + [0, 0.5];  % 调整控制点位置

        % 生成贝塞尔曲线
        t = linspace(0, 1, 40);
        t = t';

        curve = (1 - t).^2 .* repmat(startNode,length(t),1) + 2 * (1 - t) .* t .* repmat(controlPoint,length(t),1) + t.^2 .* repmat(endNode,length(t),1);
        % 绘制曲线

        plot(curve(:, 1), curve(:, 2), '--','Color',COLOR_spon, 'LineWidth', (Degree_REGION5_Total_norm(i,learnidx)+0.01)*3  );
      

    end
 

    hold on


    temp2 = temp(:,:,learnidx) - diag(diag(temp(:,:,learnidx)));
    nodeIds_chosen = [8:11]; % 节点 ID
    nodeIds_unchosen = nodeIds(~ismember(nodeIds, nodeIds_chosen)); % 节点 ID
    temp2(nodeIds_unchosen,:) = 0;
    temp2(:,nodeIds_unchosen) = 0;
%     temp2(  temp2<1.03*mean(temp,"all")  )         = 0;%1.08*mean(temp,"all")
%     temp2(  temp2<1*mean(temp,"all")  )         = 0;
    correlation_matrix_T_Mice_Learn_avg_GRAPH(:,:,learnidx)  = (temp2);

    [source, target, weight] = find(triu( temp2 ));

    edgeList = [source, target, weight];

    G = graph( source, target, weight, Cortex_13regions_new_cell );
   
  
   
    nodesize    = floor((nodeDegrees)*20+1);
   
   
   
    % 绘制曲线边
    for i = 1:size(source, 1)
        % 获取边的起点和终点

        if(  ismember( source(i),nodeIds_chosen )  )
           if(ismember( target(i),nodeIds_chosen ))
             startNode = nodes(source(i), :);
             endNode   = nodes(target(i), :);
           end
         end
        % 计算控制点（用于贝塞尔曲线）
        controlPoint = (startNode + endNode) / 2 + [0, 0.5];  % 调整控制点位置

        % 生成贝塞尔曲线
        t = linspace(0, 1, 40);
        t = t';

        curve = (1 - t).^2 .* repmat(startNode,length(t),1) + 2 * (1 - t) .* t .* repmat(controlPoint,length(t),1) + t.^2 .* repmat(endNode,length(t),1);
        % 绘制曲线

        plot(curve(:, 1), curve(:, 2), '--','Color',COLOR_spon, 'LineWidth', (Degree_REGION5_Total_norm(i+length(weight1),learnidx)+0.01)   );
      

    end

    hold on

  
    plot([3.05 4.45], [5.5 4.75],'--','Color',COLOR_spon , 'LineWidth', (Degree_REGION5_Total_norm(size(source, 1)+length(weight1)+1,learnidx)+0.01)   );
    hold on
    plot([2.55 1.25], [5.25 3], '--','Color',COLOR_spon, 'LineWidth', (Degree_REGION5_Total_norm(size(source, 1)+length(weight1)+2,learnidx)+0.01)   );
    hold on
    plot([2.55 1.25], [5.25 2], '--','Color',COLOR_spon, 'LineWidth', (Degree_REGION5_Total_norm(size(source, 1)+length(weight1)+3,learnidx)+0.01)   );
    hold on
    plot([2.55 2.25], [5.25 1.45], '--','Color',COLOR_spon, 'LineWidth', (Degree_REGION5_Total_norm(size(source, 1)+length(weight1)+4,learnidx)+0.01)   );
    hold on
   
    
    plot([4.45 3.75], [2.25 0.85],'--','Color', COLOR_spon, 'LineWidth', (Degree_REGION5_Total_norm(size(source, 1)+length(weight1)+5,learnidx)+0.01)   );
    hold on
   
    plot([2.95 1.25], [3.5 3], '--','Color',COLOR_spon, 'LineWidth', (Degree_REGION5_Total_norm(size(source, 1)+length(weight1)+6,learnidx)+0.01)   );
    hold on
    plot([2.95 1.25], [3.5 2],'--', 'Color',COLOR_spon, 'LineWidth', (Degree_REGION5_Total_norm(size(source, 1)+length(weight1)+7,learnidx)+0.01)   );
    hold on
   

    plot([0.75 0.75], [2.75 2.25], '--','Color',COLOR_spon, 'LineWidth', (Degree_REGION5_Total_norm(size(source, 1)+length(weight1)+8,learnidx)+0.01)   );
    hold on
    plot([1.25 2.25], [3 1.45], '--','Color',COLOR_spon, 'LineWidth', (Degree_REGION5_Total_norm(size(source, 1)+length(weight1)+9,learnidx)+0.01)   );
    hold on
   

    plot([1.25 2.25], [2 1.45],'--','Color', COLOR_spon, 'LineWidth', (Degree_REGION5_Total_norm(size(source, 1)+length(weight1)+10,learnidx)+0.01)   );
    hold on




  rectangle('Position', [2.8-0.75, 5.5-0.25, 1, 0.5], ...
           'FaceColor', 'none', ...
           'EdgeColor', color_inuse_GROUP(1,:), ...
           'LineWidth', 2, ...
           'Curvature', [0.2, 0.2]); % 圆角矩形

 hold on;
 

    rectangle('Position', [3.5-0.25, 2.5-0.25, 3.0, 2.5], ...
           'FaceColor', 'none', ...
           'EdgeColor', color_inuse_GROUP(4,:), ...
           'LineWidth', 2, ...
           'Curvature', [0.2, 0.2]); % 圆角矩形
     hold on;
 

    rectangle('Position', [1.0-0.75, 2.0-0.25, 1, 0.5], ...
           'FaceColor', 'none', ...
           'EdgeColor', color_inuse_GROUP(2,:), ...
           'LineWidth', 2, ...
           'Curvature', [0.2, 0.2]); % 圆角矩形
     hold on;
 

    rectangle('Position', [1.0-0.75, 3.0-0.25, 1, 0.5], ...
           'FaceColor', 'none', ...
           'EdgeColor', color_inuse_GROUP(5,:), ...
           'LineWidth', 2, ...
           'Curvature', [0.2, 0.2]); % 圆角矩形
     hold on;
 

    rectangle('Position', [1.0-0.25, 0.5-0.25, 3.0, 1.2], ...
           'FaceColor', 'none', ...
           'EdgeColor', color_inuse_GROUP(3,:), ...
           'LineWidth', 2, ...
           'Curvature', [0.2, 0.2]); % 圆角矩形
     hold on;
 



    
    scatter(x, y, nodesize*15, [color_inuse_GROUP(1,:);
        color_inuse_GROUP(4,:);
        color_inuse_GROUP(4,:);
        color_inuse_GROUP(4,:);
        color_inuse_GROUP(4,:);
        color_inuse_GROUP(4,:);
        color_inuse_GROUP(4,:);
        color_inuse_GROUP(3,:);
        color_inuse_GROUP(3,:);
        color_inuse_GROUP(3,:);
        color_inuse_GROUP(3,:);
        color_inuse_GROUP(2,:);
        color_inuse_GROUP(5,:)
        ] , 'filled');

     %text(x, y, G.Nodes.Name, 'VerticalAlignment','bottom', 'HorizontalAlignment','right');
    xlim([0 6])
    ylim([0 7])
    set(gca, 'XColor', 'none');
    set(gca, 'XTickLabel', []);
    set(gca, 'YColor', 'none');
    set(gca, 'YTickLabel', []);

   



correlation_matrix_T_Mice_Learn_avg_GRAPH = []
Degree_REGION5_learn=[]
temp_rescale= rescale(temp);

weight2 = []
weight1 = []

for learnidx = 1:6

 
    temp2 = temp(:,:,learnidx) - diag(diag(temp(:,:,learnidx)));
    nodeIds_chosen   = [2:8]; 
    nodeIds_unchosen = nodeIds(~ismember(nodeIds, nodeIds_chosen)); % 节点 ID
    temp2(nodeIds_unchosen,:) = 0;
%     temp2(  temp2<0.5*mean(temp,"all")  )         = 0;%1.08*mean(temp,"all")%0.046
    [source, target, weight] = find(triu( temp2 ));   
    weight1(:,learnidx)     = (weight);


    temp2 = temp(:,:,learnidx) - diag(diag(temp(:,:,learnidx)));
    nodeIds_chosen = [9:11];
    nodeIds_unchosen = nodeIds(~ismember(nodeIds, nodeIds_chosen)); % 节点 ID
    temp2(nodeIds_unchosen,:) = 0;   
%     temp2(  temp2<0.5*mean(temp,"all")  )         = 0;
    [source, target, weight] = find(triu( temp2 ));  
    weight2(:,learnidx)     = (weight);
   

    Degree_MOs_RSPd                    = mean( temp(1,13,learnidx) , 'all');
    Degree_MOs_RSPagl                  = mean( temp(1,12,learnidx) , 'all');
    Degree_MOs_SSP                     = mean( temp(1,2:8,learnidx) , 'all');
    Degree_MOs_VIS                     = mean( temp(1,9:11,learnidx) , 'all');
   
    Degree_SSP_RSPd                    = mean( temp(2:8,13,learnidx) , 'all');
    Degree_SSP_RSPagl                  = mean( temp(2:8,12,learnidx) , 'all');
    Degree_SSP_VIS                     = mean( temp(2:8,9:11,learnidx) , 'all');

    Degree_RSPd_RSPagl                 = mean( temp(13,12,learnidx) , 'all');
    Degree_RSPd_VIS                    = mean( temp(13,9:11,learnidx) , 'all');

    Degree_RSPagl_VIS                  = mean( temp(12,9:11,learnidx) , 'all');

    Degree_REGION5     = ([Degree_MOs_SSP Degree_MOs_RSPd Degree_MOs_RSPagl  Degree_MOs_VIS ...
        Degree_SSP_VIS Degree_SSP_RSPd Degree_SSP_RSPagl  ...
        Degree_RSPd_RSPagl Degree_RSPd_VIS ...
        Degree_RSPagl_VIS]);

Degree_REGION5_learn(:,learnidx) = Degree_REGION5;
   

end



Degree_REGION5_Total_norm=[]
Degree_REGION5_Total=[]
Degree_REGION5_Total2=[]
Degree_REGION5_Total = [weight1;weight2;Degree_REGION5_learn]
Degree_REGION5_Total2 = Degree_REGION5_Total.^5;
TEMP = Degree_REGION5_Total2(:,2:6);
Degree_REGION5_Total_norm(:,1) = 4*rescale(Degree_REGION5_Total2(:,1));
Degree_REGION5_Total_norm(:,2:6) = 4*(TEMP - min(TEMP(:))) / (max(TEMP(:)) - min(TEMP(:)));

% Degree_REGION5_Total_norm(Degree_REGION5_Total_norm<=0.3) = nan;

% close all


% COLOR_spon = [150 255 255]/255
% COLOR_learn = [150 50 255]/255


close all

for learnidx = 1

    nodeDegrees = (diag(temp_rescale(:,:,learnidx)));

    y = [ 5.5  4.5  4.0  3.5  4.0  3.0  3.0  2.5   1.0  1.0  0.5  2.0  3.0]; % 节点的 x 坐标
    x = [ 2.8  3.5  4.5  4.0  5.5  5.5  4.0  3.5   2.8  1.5  2.5  1.0  1.0]; % 节点的 y 坐标


 
    temp2 = temp(:,:,learnidx) - diag(diag(temp(:,:,learnidx)));
    nodeIds_chosen = [2:8]; % 节点 ID
    nodeIds_unchosen = nodeIds(~ismember(nodeIds, nodeIds_chosen)); % 节点 ID
    temp2(nodeIds_unchosen,:) = 0;
%     temp2(  temp2<1.03*mean(temp,"all")  )         = 0;%1.08*mean(temp,"all")
%     temp2(  temp2<1*mean(temp,"all")  )         = 0;
    correlation_matrix_T_Mice_Learn_avg_GRAPH(:,:,learnidx)  = (temp2);

    [source, target, weight] = find(triu( temp2 ));

    edgeList = [source, target, weight];

    G = graph( source, target, weight, Cortex_13regions_new_cell );
   
 
   
    nodesize    = floor((nodeDegrees)*20+1);
   
   
    handle3=figure(3000+learnidx);
    handle3.Position = [100+learnidx*250 100 500 600];
    hold on;

    nodes = [x;y]'
    % 绘制曲线边
    for i = 1:size(source, 1)
        % 获取边的起点和终点

        if(  ismember( source(i),nodeIds_chosen )  )
           if(ismember( target(i),nodeIds_chosen ))
             startNode = nodes(source(i), :);
             endNode   = nodes(target(i), :);
           end
         end
        % 计算控制点（用于贝塞尔曲线）
        controlPoint = (startNode + endNode) / 2 + [0, 0.5];  % 调整控制点位置

        % 生成贝塞尔曲线
        t = linspace(0, 1, 20);
        t = t';

        curve = (1 - t).^2 .* repmat(startNode,length(t),1) + 2 * (1 - t) .* t .* repmat(controlPoint,length(t),1) + t.^2 .* repmat(endNode,length(t),1);
        % 绘制曲线

        plot(curve(:, 1), curve(:, 2), 'Color',COLOR_spon, 'LineWidth', (Degree_REGION5_Total_norm(i,learnidx)+0.01)   );
      

    end
 

    hold on


    temp2 = temp(:,:,learnidx) - diag(diag(temp(:,:,learnidx)));
    nodeIds_chosen = [9:11]; % 节点 ID
    nodeIds_unchosen = nodeIds(~ismember(nodeIds, nodeIds_chosen)); % 节点 ID
    temp2(nodeIds_unchosen,:) = 0;

%     temp2(  temp2<1.03*mean(temp,"all")  )      = 0;%1.08*mean(temp,"all")
%     temp2(  temp2<1*mean(temp,"all")  )         = 0;

    correlation_matrix_T_Mice_Learn_avg_GRAPH(:,:,learnidx)  = (temp2);

    [source, target, weight] = find(triu( temp2 ));

    edgeList = [source, target, weight];

    G = graph( source, target, weight, Cortex_13regions_new_cell );
   
  
   
    nodesize    = floor((nodeDegrees)*20+1);
   
   
   
    % 绘制曲线边
    for i = 1:size(source, 1)
        % 获取边的起点和终点

        if(  ismember( source(i),nodeIds_chosen )  )
           if(ismember( target(i),nodeIds_chosen ))
             startNode = nodes(source(i), :);
             endNode   = nodes(target(i), :);
           end
         end
        % 计算控制点（用于贝塞尔曲线）
        controlPoint = (startNode + endNode) / 2 + [0, 0.5];  % 调整控制点位置

        % 生成贝塞尔曲线
        t = linspace(0, 1, 20);
        t = t';

        curve = (1 - t).^2 .* repmat(startNode,length(t),1) + 2 * (1 - t) .* t .* repmat(controlPoint,length(t),1) + t.^2 .* repmat(endNode,length(t),1);
        % 绘制曲线

        plot(curve(:, 1), curve(:, 2), 'Color',COLOR_spon, 'LineWidth', (Degree_REGION5_Total_norm(i+length(weight1),learnidx)+0.01)   );
      

    end

    hold on

  
    plot([3.05 4.45], [5.5 4.75],'-','Color',COLOR_spon , 'LineWidth', (Degree_REGION5_Total_norm(size(source, 1)+length(weight1)+1,learnidx)+0.01)   );
    hold on
    plot([2.55 1.25], [5.25 3], '--','Color',COLOR_spon, 'LineWidth', (Degree_REGION5_Total_norm(size(source, 1)+length(weight1)+2,learnidx)+0.01)   );
    hold on
    plot([2.55 1.25], [5.25 2], '-','Color',COLOR_spon, 'LineWidth', (Degree_REGION5_Total_norm(size(source, 1)+length(weight1)+3,learnidx)+0.01)   );
    hold on
    plot([2.55 2.25], [5.25 1.45], '-','Color',COLOR_spon, 'LineWidth', (Degree_REGION5_Total_norm(size(source, 1)+length(weight1)+4,learnidx)+0.01)   );
    hold on
   
    
    plot([4.45 3.75], [2.25 0.85],'-','Color', COLOR_spon, 'LineWidth', (Degree_REGION5_Total_norm(size(source, 1)+length(weight1)+5,learnidx)+0.01)   );
    hold on
   
    plot([2.95 1.25], [3.5 3], '-','Color',COLOR_spon, 'LineWidth', (Degree_REGION5_Total_norm(size(source, 1)+length(weight1)+6,learnidx)+0.01)   );
    hold on
    plot([2.95 1.25], [3.5 2],'-', 'Color',COLOR_spon, 'LineWidth', (Degree_REGION5_Total_norm(size(source, 1)+length(weight1)+7,learnidx)+0.01)   );
    hold on
   

    plot([0.75 0.75], [2.75 2.25], '-','Color',COLOR_spon, 'LineWidth', (Degree_REGION5_Total_norm(size(source, 1)+length(weight1)+8,learnidx)+0.01)   );
    hold on
    plot([1.25 2.25], [3 1.45], '-','Color',COLOR_spon, 'LineWidth', (Degree_REGION5_Total_norm(size(source, 1)+length(weight1)+9,learnidx)+0.01)   );
    hold on
   

    plot([1.25 2.25], [2 1.45],'-','Color', COLOR_spon, 'LineWidth', (Degree_REGION5_Total_norm(size(source, 1)+length(weight1)+10,learnidx)+0.01)   );
    hold on




  rectangle('Position', [2.8-0.75, 5.5-0.25, 1, 0.5], ...
           'FaceColor', 'none', ...
           'EdgeColor', color_inuse_GROUP(1,:), ...
           'LineWidth', 2, ...
           'Curvature', [0.2, 0.2]); % 圆角矩形

 hold on;
 

    rectangle('Position', [3.2-0.25, 2.5-0.25, 3.0, 2.5], ...
           'FaceColor', 'none', ...
           'EdgeColor', color_inuse_GROUP(4,:), ...
           'LineWidth', 2, ...
           'Curvature', [0.2, 0.2]); % 圆角矩形
     hold on;
 

    rectangle('Position', [1.0-0.75, 2.0-0.25, 1, 0.5], ...
           'FaceColor', 'none', ...
           'EdgeColor', color_inuse_GROUP(2,:), ...
           'LineWidth', 2, ...
           'Curvature', [0.2, 0.2]); % 圆角矩形
     hold on;
 

    rectangle('Position', [1.0-0.75, 3.0-0.25, 1, 0.5], ...
           'FaceColor', 'none', ...
           'EdgeColor', color_inuse_GROUP(5,:), ...
           'LineWidth', 2, ...
           'Curvature', [0.2, 0.2]); % 圆角矩形
     hold on;
 

    rectangle('Position', [1.0-0.25, 0.5-0.25, 3.0, 1.2], ...
           'FaceColor', 'none', ...
           'EdgeColor', color_inuse_GROUP(3,:), ...
           'LineWidth', 2, ...
           'Curvature', [0.2, 0.2]); % 圆角矩形
     hold on;



    
    scatter(x, y, nodesize*15, [color_inuse_GROUP(1,:);
        color_inuse_GROUP(4,:);
        color_inuse_GROUP(4,:);
        color_inuse_GROUP(4,:);
        color_inuse_GROUP(4,:);
        color_inuse_GROUP(4,:);
        color_inuse_GROUP(4,:);
        color_inuse_GROUP(4,:);
        color_inuse_GROUP(3,:);
        color_inuse_GROUP(3,:);
        color_inuse_GROUP(3,:);
        color_inuse_GROUP(2,:);
        color_inuse_GROUP(5,:)
        ] , 'filled');

     %text(x, y, G.Nodes.Name, 'VerticalAlignment','bottom', 'HorizontalAlignment','right');
    xlim([0 6])
    ylim([0 7])
    set(gca, 'XColor', 'none');
    set(gca, 'XTickLabel', []);
    set(gca, 'YColor', 'none');
    set(gca, 'YTickLabel', []);

   



end




for learnidx = 2:6

    nodeDegrees = (diag(temp_rescale(:,:,learnidx)));

    y = [ 5.5  4.5  4.0  3.5  4.0  3.0  3.0  2.5   1.0  1.0  0.5  2.0  3.0]; % 节点的 x 坐标
    x = [ 2.8  3.5  4.5  4.0  5.5  5.5  4.0  3.5   2.8  1.5  2.5  1.0  1.0]; % 节点的 y 坐标


 


    temp2 = temp(:,:,learnidx) - diag(diag(temp(:,:,learnidx)));
    nodeIds_chosen = [2:8]; % 节点 ID
    nodeIds_unchosen = nodeIds(~ismember(nodeIds, nodeIds_chosen)); % 节点 ID
    temp2(nodeIds_unchosen,:) = 0;
%     temp2(  temp2<1.03*mean(temp,"all")  )         = 0;%1.08*mean(temp,"all")
%     temp2(  temp2<1*mean(temp,"all")  )         = 0;
    correlation_matrix_T_Mice_Learn_avg_GRAPH(:,:,learnidx)  = (temp2);

    [source, target, weight] = find(triu( temp2 ));

    edgeList = [source, target, weight];

    G = graph( source, target, weight, Cortex_13regions_new_cell );
   
 
   
    nodesize    = floor((nodeDegrees)*20+1);
   
   
    handle3=figure(3000+learnidx);
    handle3.Position = [100+learnidx*250 100 500 600];
    hold on;

    nodes = [x;y]'
    % 绘制曲线边
    for i = 1:size(source, 1)
        % 获取边的起点和终点

        if(  ismember( source(i),nodeIds_chosen )  )
           if(ismember( target(i),nodeIds_chosen ))
             startNode = nodes(source(i), :);
             endNode   = nodes(target(i), :);
           end
         end
        % 计算控制点（用于贝塞尔曲线）
        controlPoint = (startNode + endNode) / 2 + [0, 0.5];  % 调整控制点位置

        % 生成贝塞尔曲线
        t = linspace(0, 1, 100);
        t = t';

        curve = (1 - t).^2 .* repmat(startNode,length(t),1) + 2 * (1 - t) .* t .* repmat(controlPoint,length(t),1) + t.^2 .* repmat(endNode,length(t),1);
        % 绘制曲线

        plot(curve(:, 1), curve(:, 2), 'Color',COLOR_learn, 'LineWidth', (Degree_REGION5_Total_norm(i,learnidx)+0.01)   );
      

    end
 

    hold on


    temp2 = temp(:,:,learnidx) - diag(diag(temp(:,:,learnidx)));
    nodeIds_chosen = [9:11]; % 节点 ID
    nodeIds_unchosen = nodeIds(~ismember(nodeIds, nodeIds_chosen)); % 节点 ID
    temp2(nodeIds_unchosen,:) = 0;
%     temp2(  temp2<1.03*mean(temp,"all")  )         = 0;%1.08*mean(temp,"all")
%     temp2(  temp2<1*mean(temp,"all")  )         = 0;
    correlation_matrix_T_Mice_Learn_avg_GRAPH(:,:,learnidx)  = (temp2);

    [source, target, weight] = find(triu( temp2 ));

    edgeList = [source, target, weight];

    G = graph( source, target, weight, Cortex_13regions_new_cell );
   
  
   
    nodesize    = floor((nodeDegrees)*20+1);
   
   
   
    % 绘制曲线边
    for i = 1:size(source, 1)
        % 获取边的起点和终点

        if(  ismember( source(i),nodeIds_chosen )  )
           if(ismember( target(i),nodeIds_chosen ))
             startNode = nodes(source(i), :);
             endNode   = nodes(target(i), :);
           end
         end
        % 计算控制点（用于贝塞尔曲线）
        controlPoint = (startNode + endNode) / 2 + [0, 0.5];  % 调整控制点位置

        % 生成贝塞尔曲线
        t = linspace(0, 1, 100);
        t = t';

        curve = (1 - t).^2 .* repmat(startNode,length(t),1) + 2 * (1 - t) .* t .* repmat(controlPoint,length(t),1) + t.^2 .* repmat(endNode,length(t),1);
        % 绘制曲线

        plot(curve(:, 1), curve(:, 2), 'Color',COLOR_learn, 'LineWidth', (Degree_REGION5_Total_norm(i+length(weight1),learnidx)+0.01)   );
      

    end

    hold on

  
    plot([3.05 4.45], [5.5 4.75],'-','Color',COLOR_spon , 'LineWidth', (Degree_REGION5_Total_norm(size(source, 1)+length(weight1)+1,learnidx)+0.01)   );
    hold on
    plot([2.55 1.25], [5.25 3], '-','Color',COLOR_spon, 'LineWidth', (Degree_REGION5_Total_norm(size(source, 1)+length(weight1)+2,learnidx)+0.01)   );
    hold on
    plot([2.55 1.25], [5.25 2], '-','Color',COLOR_spon, 'LineWidth', (Degree_REGION5_Total_norm(size(source, 1)+length(weight1)+3,learnidx)+0.01)   );
    hold on
    plot([2.55 2.25], [5.25 1.45], '-','Color',COLOR_spon, 'LineWidth', (Degree_REGION5_Total_norm(size(source, 1)+length(weight1)+4,learnidx)+0.01)   );
    hold on
   
    
    plot([4.45 3.75], [2.25 0.85],'-','Color', COLOR_spon, 'LineWidth', (Degree_REGION5_Total_norm(size(source, 1)+length(weight1)+5,learnidx)+0.01)   );
    hold on
   
    plot([2.95 1.25], [3.5 3], '-','Color',COLOR_spon, 'LineWidth', (Degree_REGION5_Total_norm(size(source, 1)+length(weight1)+6,learnidx)+0.01)   );
    hold on
    plot([2.95 1.25], [3.5 2],'-', 'Color',COLOR_spon, 'LineWidth', (Degree_REGION5_Total_norm(size(source, 1)+length(weight1)+7,learnidx)+0.01)   );
    hold on
   

    plot([0.75 0.75], [2.75 2.25], '-','Color',COLOR_spon, 'LineWidth', (Degree_REGION5_Total_norm(size(source, 1)+length(weight1)+8,learnidx)+0.01)   );
    hold on
    plot([1.25 2.25], [3 1.45], '-','Color',COLOR_spon, 'LineWidth', (Degree_REGION5_Total_norm(size(source, 1)+length(weight1)+9,learnidx)+0.01)   );
    hold on
   

    plot([1.25 2.25], [2 1.45],'-','Color', COLOR_spon, 'LineWidth', (Degree_REGION5_Total_norm(size(source, 1)+length(weight1)+10,learnidx)+0.01)   );
    hold on




  rectangle('Position', [2.8-0.75, 5.5-0.25, 1, 0.5], ...
           'FaceColor', 'none', ...
           'EdgeColor', color_inuse_GROUP(1,:), ...
           'LineWidth', 2, ...
           'Curvature', [0.2, 0.2]); % 圆角矩形

 hold on;
 

    rectangle('Position', [3.2-0.25, 2.5-0.25, 3.0, 2.5], ...
           'FaceColor', 'none', ...
           'EdgeColor', color_inuse_GROUP(4,:), ...
           'LineWidth', 2, ...
           'Curvature', [0.2, 0.2]); % 圆角矩形
     hold on;
 

    rectangle('Position', [1.0-0.75, 2.0-0.25, 1, 0.5], ...
           'FaceColor', 'none', ...
           'EdgeColor', color_inuse_GROUP(2,:), ...
           'LineWidth', 2, ...
           'Curvature', [0.2, 0.2]); % 圆角矩形
     hold on;
 

    rectangle('Position', [1.0-0.75, 3.0-0.25, 1, 0.5], ...
           'FaceColor', 'none', ...
           'EdgeColor', color_inuse_GROUP(5,:), ...
           'LineWidth', 2, ...
           'Curvature', [0.2, 0.2]); % 圆角矩形
     hold on;
 

    rectangle('Position', [1.0-0.25, 0.5-0.25, 3.0, 1.2], ...
           'FaceColor', 'none', ...
           'EdgeColor', color_inuse_GROUP(3,:), ...
           'LineWidth', 2, ...
           'Curvature', [0.2, 0.2]); % 圆角矩形
     hold on;



%     temp2 = temp(:,:,learnidx) - diag(diag(temp(:,:,learnidx)));
%     
%     nodeDegrees = (diag(temp_rescale(:,:,learnidx)));
% 
% 
%     temp2(  temp2<0.1  )         = 0;%1.08*mean(temp,"all")%0.046
% 
%     correlation_matrix_T_Mice_Learn_avg_GRAPH(:,:,learnidx)  = (temp2);
% 
%     [source, target, weight] = find(triu( correlation_matrix_T_Mice_Learn_avg_GRAPH(:,:,learnidx) ));
% 
%     edgeList = [source, target, weight];
% 
%     G = graph( source, target, weight, Cortex_13regions_new_cell );
% 
% 
%     y = [ 5.5  4.5  4.0  3.5  4.0  3.0  3.0  2.5   1.0  1.0  0.5  2.0  3.0]; % 节点的 x 坐标
%     x = [ 2.8  3.5  4.5  4.0  5.5  5.5  4.0  3.5   2.8  1.5  2.5  1.0  1.0]; % 节点的 y 坐标
% 
% 
%     weight2 =rescale(weight);
%    
%     nodesize    = floor((nodeDegrees)*10+1);
%     % 可视化图
%    
%     handle3=figure(1100+learnidx);
%     handle3.Position = [100+learnidx*250 100 500 600];
%     hold on;
% 
%     nodes = [x;y]'
%     % 绘制曲线边
%     for i = 1:size(source, 1)
%         % 获取边的起点和终点
%         startNode = nodes(source(i), :);
%         endNode   = nodes(target(i), :);
% 
%         % 计算控制点（用于贝塞尔曲线）
%         controlPoint = (startNode + endNode) / 2 + [0, 0.5];  % 调整控制点位置
% 
%         % 生成贝塞尔曲线
%         t = linspace(0, 1, 20);
%         t = t';
% 
%         curve = (1 - t).^2 .* repmat(startNode,length(t),1) + 2 * (1 - t) .* t .* repmat(controlPoint,length(t),1) + t.^2 .* repmat(endNode,length(t),1);
%         % 绘制曲线
% 
% %         plot(curve(:, 1), curve(:, 2), 'k.', 'LineWidth',  0.1 );%(weight2(i)+0.01)*1
%  plot(curve(:, 1), curve(:, 2), 'k', 'LineWidth',  (weight2(i)+0.01)*2 );%(weight2(i)+0.01)*1
%     end
% 
%     hold on;
%  

% %  rectangle('Position', [2.8-0.25, 5.5-0.25, 0.5, 0.5], ...
% %            'FaceColor', 'none', ...
% %            'EdgeColor', 'red', ...
% %            'LineWidth', 1, ...
% %            'Curvature', [0.2, 0.2]); % 圆角矩形
% % 
% %  hold on;
% %  
% % 
% %     rectangle('Position', [3.2-0.25, 2.5-0.25, 3.0, 2.5], ...
% %            'FaceColor', 'none', ...
% %            'EdgeColor', 'red', ...
% %            'LineWidth', 1, ...
% %            'Curvature', [0.2, 0.2]); % 圆角矩形
% %      hold on;
% %  
% % 
% %     rectangle('Position', [1.0-0.25, 2.0-0.25, 0.5, 0.5], ...
% %            'FaceColor', 'none', ...
% %            'EdgeColor', 'red', ...
% %            'LineWidth', 1, ...
% %            'Curvature', [0.2, 0.2]); % 圆角矩形
% %      hold on;
% %  
% % 
% %     rectangle('Position', [1.0-0.25, 3.0-0.25, 0.5, 0.5], ...
% %            'FaceColor', 'none', ...
% %            'EdgeColor', 'red', ...
% %            'LineWidth', 1, ...
% %            'Curvature', [0.2, 0.2]); % 圆角矩形
% %      hold on;
% %  
% % 
% %     rectangle('Position', [1.0-0.25, 0.5-0.25, 2.5, 1.2], ...
% %            'FaceColor', 'none', ...
% %            'EdgeColor', 'red', ...
% %            'LineWidth', 1, ...
% %            'Curvature', [0.2, 0.2]); % 圆角矩形
% %      hold on;
% %  



    
    scatter(x, y, nodesize*15, [color_inuse_GROUP(1,:);
        color_inuse_GROUP(4,:);
        color_inuse_GROUP(4,:);
        color_inuse_GROUP(4,:);
        color_inuse_GROUP(4,:);
        color_inuse_GROUP(4,:);
        color_inuse_GROUP(4,:);
        color_inuse_GROUP(4,:);
        color_inuse_GROUP(3,:);
        color_inuse_GROUP(3,:);
        color_inuse_GROUP(3,:);
        color_inuse_GROUP(2,:);
        color_inuse_GROUP(5,:)
        ] , 'filled');

     %text(x, y, G.Nodes.Name, 'VerticalAlignment','bottom', 'HorizontalAlignment','right');
    xlim([0 6])
    ylim([0 7])
    set(gca, 'XColor', 'none');
    set(gca, 'XTickLabel', []);
    set(gca, 'YColor', 'none');
    set(gca, 'YTickLabel', []);

   



end





nodeIds        = [1; 2;3;4;5;6;7;8;9;10;11;12;13]; % 节点 ID
nodeIds_unchosen = [2;3;4;5;6;7;13]; % 节点 ID
nodeIds_chosen = [1;8;9;10;11;12]; % 节点 ID
nodeLabels     = Cortex_13regions_new_cell(nodeIds); % 节点标签
nodeCategories = {'Group1'; 'Group2'; 'Group2'; 'Group2'; 'Group2'; 'Group2'; 'Group2'; 'Group2';'Group3';'Group3';'Group3';'Group4';'Group5'}; % 节点类别
nodeCategories_chosen = nodeCategories(nodeIds)

color_inuse_GROUP                  = []
color_inuse_GROUP                  = color_inuse([19 3 2 16 24],:);
color_inuse_GROUP(6,:)             = [0.5 0.5 0.5];


temp                               = double(DegreeCentrality_learn_Reg_MiceAVG_13Regions);%>mean(DegreeCentrality_learn_Reg_MiceAVG_13Regions(:,:,learnidx),'all'));
correlation_matrix_T_Mice_Learn_avg_GRAPH = []

temp_rescale= rescale(temp);

close all

for learnidx = 1:6


    temp2 = temp(:,:,learnidx) - diag(diag(temp(:,:,learnidx)));
    temp2(nodeIds_unchosen,:) = 0;


    nodeDegrees = (diag(temp_rescale(:,:,learnidx)));


    temp2(  temp2<1.08*mean(temp,"all")  )         = 0;%1.08*mean(temp,"all")

    correlation_matrix_T_Mice_Learn_avg_GRAPH(:,:,learnidx)  = (temp2);

    [source, target, weight] = find(triu( correlation_matrix_T_Mice_Learn_avg_GRAPH(:,:,learnidx) ));

    edgeList = [source, target, weight];

    G = graph( source, target, weight, Cortex_13regions_new_cell );

    y = [ 5  5  4  3.5  4  2.6  2.6  1.6  1.2  1  1  2  2]; % 节点的 x 坐标
    x = [ 2  4  1  3    5  2    4    3    3    2  4  1  5]; % 节点的 y 坐标

    y = [ 6  5  4  3.5  4  2    2.6  1.6  1.2  1  1  2  2]; % 节点的 x 坐标
    x = [ 3  3  1  3    5  3    3    3    3    2  4  1  5]; % 节点的 y 坐标

    y = [ 6  5    3.5  4  4  2.5  3.5  2   1.2  1  1  2  2]; % 节点的 x 坐标
    x = [ 3  2.8  2    1  5  3.3  4    2.8   3.2  2  4  1  5]; % 节点的 y 坐标


    weight2     = rescale(weight);
   
    nodesize    = floor((nodeDegrees)*20+1);
    % 可视化图
   
    handle3=figure(1000+learnidx);
    handle3.Position = [100+learnidx*250 100 500 600];
    hold on;

    nodes = [x;y]'
    % 绘制曲线边
    for i = 1:size(source, 1)
        % 获取边的起点和终点

        if(  ismember( source(i),nodeIds_chosen )  )
           if(ismember( target(i),nodeIds_chosen ))
             startNode = nodes(source(i), :);
             endNode   = nodes(target(i), :);
           end
         end
        % 计算控制点（用于贝塞尔曲线）
        controlPoint = (startNode + endNode) / 2 + [0, 0.5];  % 调整控制点位置

        % 生成贝塞尔曲线
        t = linspace(0, 1, 100);
        t = t';

        curve = (1 - t).^2 .* repmat(startNode,length(t),1) + 2 * (1 - t) .* t .* repmat(controlPoint,length(t),1) + t.^2 .* repmat(endNode,length(t),1);
        % 绘制曲线

        plot(curve(:, 1), curve(:, 2), 'k', 'LineWidth', (weight2(i)+0.01)  );
      

    end

    hold on;
    scatter(x(nodeIds_chosen), y(nodeIds_chosen), nodesize(nodeIds_chosen)*15, [color_inuse_GROUP(1,:);
       
        color_inuse_GROUP(4,:);
        color_inuse_GROUP(3,:);
        color_inuse_GROUP(3,:);
        color_inuse_GROUP(3,:);
        color_inuse_GROUP(2,:);
    
        ] , 'filled');

%     text(x, y, G.Nodes.Name, 'VerticalAlignment','bottom', 'HorizontalAlignment','right');
    xlim([0 6])
    ylim([0 7])
    set(gca, 'XColor', 'none');
    set(gca, 'XTickLabel', []);
    set(gca, 'YColor', 'none');
    set(gca, 'YTickLabel', []);

end







temp      = double(DegreeCentrality_cue_Reg_MiceAVG_13Regions);%>mean(DegreeCentrality_learn_Reg_MiceAVG_13Regions(:,:,learnidx),'all'));
correlation_matrix_T_Mice_Learn_avg_GRAPH = []



for learnidx = 1:5


      temp2 = temp(:,:,learnidx) - diag(diag(temp(:,:,learnidx)));
    
    nodeDegrees = (diag(temp_rescale(:,:,learnidx)));


    temp2(  temp2<1*mean(temp,"all")  )         = 0;

    correlation_matrix_T_Mice_Learn_avg_GRAPH(:,:,learnidx)  = (temp2);

    [source, target, weight] = find(triu( correlation_matrix_T_Mice_Learn_avg_GRAPH(:,:,learnidx) ));

    edgeList = [source, target, weight];

    G = graph( source, target, weight, Cortex_13regions_new_cell );

    y = [ 5  5  4  3.5  4  2.6  2.6  1.6  1.2  1  1  2  2]; % 节点的 x 坐标
    x = [ 2  4  1  3    5  2    4    3    3    2  4  1  5]; % 节点的 y 坐标

    y = [ 6  5  4  3.5  4  2    2.6  1.6  1.2  1  1  2  2]; % 节点的 x 坐标
    x = [ 3  3  1  3    5  3    3    3    3    2  4  1  5]; % 节点的 y 坐标

    y = [ 6  5    3.5  4  4  2.5  3.5  1.5   1.2  1  1  2  2]; % 节点的 x 坐标
    x = [ 3  2.8  2    1  5  3.3  4    2.8   3.2  2  4  1  5]; % 节点的 y 坐标


    weight2 =rescale(weight);
   
    nodesize    = floor((nodeDegrees)*10+1);
    % 可视化图
   
    handle3=figure(3000+learnidx);
    handle3.Position = [100 100 500 600];
    hold on;

    nodes = [x;y]'
    % 绘制曲线边
    for i = 1:size(source, 1)
        % 获取边的起点和终点
        startNode = nodes(source(i), :);
        endNode   = nodes(target(i), :);

        % 计算控制点（用于贝塞尔曲线）
        controlPoint = (startNode + endNode) / 2 + [0, 0.5];  % 调整控制点位置

        % 生成贝塞尔曲线
        t = linspace(0, 1, 100);
        t = t';

        curve = (1 - t).^2 .* repmat(startNode,length(t),1) + 2 * (1 - t) .* t .* repmat(controlPoint,length(t),1) + t.^2 .* repmat(endNode,length(t),1);
        % 绘制曲线

        plot(curve(:, 1), curve(:, 2), 'k', 'LineWidth', (weight2(i)+0.01)*2);

    end

    hold on;
    scatter(x, y, nodesize*15, [color_inuse_GROUP(1,:);
        color_inuse_GROUP(4,:);
        color_inuse_GROUP(4,:);
        color_inuse_GROUP(4,:);
        color_inuse_GROUP(4,:);
        color_inuse_GROUP(4,:);
        color_inuse_GROUP(4,:);
        color_inuse_GROUP(4,:);
        color_inuse_GROUP(3,:);
        color_inuse_GROUP(3,:);
        color_inuse_GROUP(3,:);
        color_inuse_GROUP(2,:);
        color_inuse_GROUP(5,:)
        ] , 'filled');
text(x, y, G.Nodes.Name, 'VerticalAlignment','bottom', 'HorizontalAlignment','right');
    xlim([0 6])
    ylim([0 7])
    set(gca, 'XColor', 'none');
    set(gca, 'XTickLabel', []);
    set(gca, 'YColor', 'none');
    set(gca, 'YTickLabel', []);

end



%%-------------Cortex regions degree analysis------------------------------
%%=========================================================================
%%=========================================================================
Cortex_Incommon_New       = Cortex_Incommon([3:4, 1:2, 15:16,  11:12,   17:18,9:10,13:14, 19:20,    21:22, 25:26, 23:24,5:6, 7:8]);
[~, indexInArray2]        = ismember(Cortex_Incommon_New,Cortex_Incommon_Sorted);



Frame_chosen_index        = Frame_Learn_index;

close all

G_neurons_learn           = {}


for mouseidx              = 1:5


    Trace_Region_1        = []
    Trace_Region_2        = []
    Trace_Region_T        = []

    regionidx1            = [9:12];

    temp1                 = find(ismember(Mice_STRUCT_T(mouseidx).Neuron_REGION_ID, CortexID_Mouse_Sorted(mouseidx,indexInArray2(regionidx1))));

    Trace_Region_1(:,:)   = Mice_STRUCT_T(mouseidx).NEURON.Trace_unify(temp1 , logical(Frame_Inuse_Mouse{mouseidx}));

    Position_1            = Mice_STRUCT_T(mouseidx).COORDINATES(temp1,:);

    x1                    = Position_1(:,1);
    y1                    = Position_1(:,2);

    regionidx2            = [17:22];

    temp2                 = find(ismember(Mice_STRUCT_T(mouseidx).Neuron_REGION_ID, CortexID_Mouse_Sorted(mouseidx,indexInArray2(regionidx2))));

    Trace_Region_2(:,:)   = Mice_STRUCT_T(mouseidx).NEURON.Trace_unify(temp2 , logical(Frame_Inuse_Mouse{mouseidx}));

    Position_2            = Mice_STRUCT_T(mouseidx).COORDINATES(temp2,:);

    x2                    = Position_2(:,1);
    y2                    = Position_2(:,2);

    x_T                   = [x1;x2];
    y_T                   = [y1;y2];

    Trace_Region_T        = [Trace_Region_1;Trace_Region_2];


    degree_neurons        = []

    for learnidx          = 1:max(Frame_chosen_index)

        Trace_Region_T_learn       = Trace_Region_T(:,Frame_chosen_index==learnidx);

        correlation_matrix_T       = corr(Trace_Region_T_learn');
        correlation_matrix_T(isnan(correlation_matrix_T)) = 0;

        binary_matrix_T            = correlation_matrix_T ;
        binary_matrix_T(binary_matrix_T<0.9) = 0;

        binary_matrix_T_adjcent    = binary_matrix_T - diag(diag(binary_matrix_T));

        %         [source1, target1, weight1]   = find(triu( binary_matrix_T ));
        %         G_neurons_learn_temp          = graph( source, target, weight );

        [source, target, weight]   = find(triu( binary_matrix_T ));

        G_neurons_learn{learnidx}  = graph( source, target, weight );

        G_neurons_learn{learnidx}.Nodes.X = x_T; % X坐标
        G_neurons_learn{learnidx}.Nodes.Y = y_T; % Y坐标

        figure(6000+learnidx)

        for regionidx=1:length(Mice_STRUCT_T(mouseidx).REGION_ID)

            BW         = imbinarize(double(Mice_STRUCT_T(mouseidx).ATLAS2COMET ==regionidx));
            boundaries = bwboundaries(BW);

            if (~isempty(boundaries))
                for k  = 1:1

                    b  = boundaries{k};
                    b(:,2) = b(:,2) + 985 - Bregma_Position(mouseidx,1);

                    plot(b(:,2), b(:,1), 'Color' , [0 , 0 , 0]/255, 'LineWidth', 1); % Plot boundary with green line
                    %     fill(b(:,2)+1, b(:,1)+1, Mice_STRUCT_T(1).COLOR(11,:),'FaceAlpha',0.5,'LineStyle',':');
                    hold on

                end
                hold on
            end
            hold on
            regionidx
        end


        p=plot( G_neurons_learn{learnidx}, 'XData', x_T+ 985 - Bregma_Position(mouseidx,1), 'YData', y_T, 'EdgeAlpha', 0.1);%,'EdgeLabel', round(G.Edges.Weight,1)
        p.NodeColor=[ 0 0 0.8  ];

        p.EdgeColor=[255 100 100]/255;
        p.LineWidth=0.1;


        selfLoopIdx = find(diag(binary_matrix_T)); % 找到自环边的索引
        if ~isempty(selfLoopIdx)
            highlight(p, selfLoopIdx, selfLoopIdx, 'LineStyle', 'none');
        end
        axis equal


        %         title([Cortex_13regions_new_cell{regionidx} '-M: ' num2str(mouseidx) '-Learn: ' num2str(learnidx)])

        hold on

        %         pause

    end

    pause

end




Frame_chosen_index        = Frame_Learn_index;

close all

G_neurons_learn           = {}


for mouseidx              = 1%1:5


    Trace_Region_1        = []
    Trace_Region_2        = []
    Trace_Region_T        = []

    regionidx1            = [1];

    temp1                 = find(ismember(Mice_STRUCT_T(mouseidx).Neuron_REGION_ID, CortexID_Mouse_Sorted(mouseidx,indexInArray2(regionidx1))));

    Trace_Region_1(:,:)   = Mice_STRUCT_T(mouseidx).NEURON.Trace_unify(temp1 , logical(Frame_Inuse_Mouse{mouseidx}));

    Position_1            = Mice_STRUCT_T(mouseidx).COORDINATES(temp1,:);

    x1                    = Position_1(:,1);
    y1                    = Position_1(:,2);

    regionidx2            = [2:26];

    temp2                 = find(ismember(Mice_STRUCT_T(mouseidx).Neuron_REGION_ID, CortexID_Mouse_Sorted(mouseidx,indexInArray2(regionidx2))));

    Trace_Region_2(:,:)   = Mice_STRUCT_T(mouseidx).NEURON.Trace_unify(temp2 , logical(Frame_Inuse_Mouse{mouseidx}));

    Position_2            = Mice_STRUCT_T(mouseidx).COORDINATES(temp2,:);

    x2                    = Position_2(:,1);
    y2                    = Position_2(:,2);

    x_T                   = [x1;x2];
    y_T                   = [y1;y2];

    Trace_Region_T        = [Trace_Region_1;Trace_Region_2];


    degree_neurons        = []

    for learnidx          = 1:max(Frame_chosen_index)

        Trace_Region_T_learn       = Trace_Region_T(:,Frame_chosen_index==learnidx);

        correlation_matrix_T       = corr(Trace_Region_T_learn');
        correlation_matrix_T(isnan(correlation_matrix_T)) = 0;

        binary_matrix_T            = correlation_matrix_T ;
        binary_matrix_T(binary_matrix_T<0.9) = 0;

        binary_matrix_T_adjcent    = binary_matrix_T - diag(diag(binary_matrix_T));

        nodeDegrees = sum(binary_matrix_T_adjcent,2);

        %         [source1, target1, weight1]   = find(triu( binary_matrix_T ));
        %         G_neurons_learn_temp          = graph( source, target, weight );

        [source, target, weight]   = find(triu( binary_matrix_T ));

        G_neurons_learn{learnidx}  = graph( source, target, weight );

        G_neurons_learn{learnidx}.Nodes.X = x_T; % X坐标
        G_neurons_learn{learnidx}.Nodes.Y = y_T; % Y坐标

        figure(6000+learnidx)

        for regionidx=1:length(Mice_STRUCT_T(mouseidx).REGION_ID)

            BW         = imbinarize(double(Mice_STRUCT_T(mouseidx).ATLAS2COMET ==regionidx));
            boundaries = bwboundaries(BW);

            if (~isempty(boundaries))
                for k  = 1:1

                    b  = boundaries{k};
                    b(:,2) = b(:,2) + 985 - Bregma_Position(mouseidx,1);

                    plot(b(:,2), b(:,1), 'Color' , [0 , 0 , 0]/255, 'LineWidth', 1); % Plot boundary with green line
                    %     fill(b(:,2)+1, b(:,1)+1, Mice_STRUCT_T(1).COLOR(11,:),'FaceAlpha',0.5,'LineStyle',':');
                    hold on

                end
                hold on
            end
            hold on
            regionidx
        end


        p=plot( G_neurons_learn{learnidx}, 'XData', x_T+ 985 - Bregma_Position(mouseidx,1), 'YData', y_T, 'NodeColor', 'none', 'EdgeAlpha', 0.0);%,'EdgeLabel', round(G.Edges.Weight,1)
      
        p.EdgeColor=repmat([255 255 255]/255, numedges(G_neurons_learn{learnidx}), 1);
        p.LineWidth=0.1;


        selfLoopIdx = find(diag(binary_matrix_T)); % 找到自环边的索引
        if ~isempty(selfLoopIdx)
            highlight(p, selfLoopIdx, selfLoopIdx, 'LineStyle', 'none');
        end
        axis equal


        %         title([Cortex_13regions_new_cell{regionidx} '-M: ' num2str(mouseidx) '-Learn: ' num2str(learnidx)])

        hold on

              scatter(x_T+ 985 - Bregma_Position(mouseidx,1), y_T, 5, [0.5 0.5 0.5],  'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'none' , 'LineWidth', 0.1);
              hold on;

        targetNodes       = find(nodeDegrees>20)';
        hold on
        % 标记这些节点的所有边（红色）
        for node = targetNodes
            %         highlight(p , 'Edges', outedges( G_neurons_learn{learnidx}, node), 'EdgeColor', [178 102 255]/255);
            highlight(p , node, 'NodeColor', 'r'); % 节点设为蓝色
            hold on
            edgeIndices = outedges( G_neurons_learn{learnidx}, node);

            p.EdgeColor(edgeIndices, :) = repmat([1 0 0], length(edgeIndices), 1);
            p.EdgeAlpha = 0.2;
            hold on

        end


        %         pause

    end

    pause

end





[degree_neurons_sort index_sort] = sortrows(degree_neurons, 3);
degree_neurons_diff              = diff(degree_neurons,1,2);
degree_neurons_diff_b            = diff(degree_neurons,1,2)>0;

degree_neurons_diff_increase_index = find((sum(degree_neurons_diff,2)>=5)&(degree_neurons(:,2)>=100))
degree_neurons_diff_increase_list  = ((sum(degree_neurons_diff,2)>=5)&(degree_neurons(:,2)>=100))


close all

for learnidx=1:6

    figure(6200+learnidx)


    nodeDegrees       = degree( G_neurons_learn{learnidx} );
    DegreeIncrNodes   = degree_neurons_diff_increase_index;


    nodeAlphaValues   = ( degree_neurons_diff_increase_list );
    nodeDegrees_scale = rescale(nodeDegrees,0.1 , 1);
    subG = subgraph(G_neurons_learn{learnidx}, DegreeIncrNodes);

    x2 = G_neurons_learn{learnidx}.Nodes.X ; % X坐标
    y2 = G_neurons_learn{learnidx}.Nodes.Y ; % Y坐标


    p                 = plot( G_neurons_learn{learnidx},   'MarkerSize', 0.1,  'XData', x2+ 985 - Bregma_Position(mouseidx,1), 'YData', y2, 'EdgeColor', [178 102 255]/255 , 'EdgeAlpha', 0.0);%,'EdgeLabel', round(G.Edges.Weight,1)


    targetNodes       = find(nodeDegrees>100)';
    hold on
    % 标记这些节点的所有边（红色）
    for node = targetNodes
        highlight(p , 'Edges', outedges( G_neurons_learn{learnidx}, node), 'EdgeColor', [178 102 255]/255);
        highlight(p , node, 'NodeColor', 'r'); % 节点设为蓝色
        hold on
    end



    hold on
    scatter(x, y, 5, [0.5 0.5 0.5], 'filled','AlphaData', 0.3 );

    %     hold on
    %     scatter(x(nodeAlphaValues), y(nodeAlphaValues), 50*(nodeDegrees_scale(nodeAlphaValues)), 'r', 'filled' );



    axis equal


    %     hold on
    %
    %     x              = Position_T(:,1);
    %     y              = Position_T(:,2);
    %     for regionidx=1:length(Mice_STRUCT_T(mouseidx).REGION_ID)
    %
    %         BW         = imbinarize(double(Mice_STRUCT_T(mouseidx).ATLAS2COMET ==regionidx));
    %         boundaries = bwboundaries(BW);
    %
    %         if (~isempty(boundaries))
    %             for k  = 1:1
    %
    %                 b  = boundaries{k};
    %                 b(:,2) = b(:,2) + 985 - Bregma_Position(mouseidx,1);
    %
    %                 plot(b(:,2), b(:,1), 'Color' , [0 , 0 , 0]/255, 'LineWidth', 1); % Plot boundary with green line
    %                 %     fill(b(:,2)+1, b(:,1)+1, Mice_STRUCT_T(1).COLOR(11,:),'FaceAlpha',0.5,'LineStyle',':');
    %                 hold on
    %
    %             end
    %             hold on
    %         end
    %         hold on
    %         regionidx
    %     end
    %
    %
    %     hold off
    %
        %pause


end






close all

for learnidx=1:6

    figure(6500+learnidx)



    nodeDegrees       = degree( G_neurons_learn{learnidx} );
    highDegreeNodes   = find(nodeDegrees > 100);

    %     [s, t]            = findedge( G_neurons_learn{learnidx} ); % 获取所有边的起点和终点
    % %     edgesToShow       = ismember(s, highDegreeNodes) | ismember(t, highDegreeNodes); % 标记需要显示的边
    % %     edgeColors        = repmat(color_inuse_GROUP(learnidx,:), numedges( G_neurons_learn{learnidx}), 1); % 默认黑色
    % %     edgeColors(~edgesToShow, :) = repmat([1, 1, 1], sum(~edgesToShow), 1); % 隐藏不需要的边（设置为白色）


    nodeAlphaValues   = (nodeDegrees > 100);
    nodeDegrees_scale = rescale(nodeDegrees,0.1 , 1);
    subG = subgraph(G_neurons_learn{learnidx}, highDegreeNodes);

    x2 = G_neurons_learn{learnidx}.Nodes.X ; % X坐标
    y2 = G_neurons_learn{learnidx}.Nodes.Y ; % Y坐标



    p                 = plot( subG,   'MarkerSize', 0.1,  'XData', x2(highDegreeNodes)+ 985 - Bregma_Position(mouseidx,1), 'YData', y(highDegreeNodes), 'EdgeColor', [178 102 255]/255 , 'EdgeAlpha', 0.2);%,'EdgeLabel', round(G.Edges.Weight,1)

    %     edges = G_neurons_learn{learnidx}.Edges.EndNodes; % 获取所有边
    %     xline = [G_neurons_learn{learnidx}.Nodes.X(edges(:,1)), G_neurons_learn{learnidx}.Nodes.X(edges(:,2))]'; % 起点和终点的X坐标
    %     yline = [G_neurons_learn{learnidx}.Nodes.Y(edges(:,1)), G_neurons_learn{learnidx}.Nodes.Y(edges(:,2))]'; % 起点和终点的Y坐标
    %     line(xline, yline, 'Color', 'b'); % 绘制所有边
    hold on
    scatter(x, y, 5, [0.5 0.5 0.5], 'filled','AlphaData', 0.3 );

    hold on
    scatter(x(nodeAlphaValues), y(nodeAlphaValues), 50*(nodeDegrees_scale(nodeAlphaValues)), 'r', 'filled' );




    axis equal


    hold on

    x              = Position_T(:,1);
    y              = Position_T(:,2);
    for regionidx=1:length(Mice_STRUCT_T(mouseidx).REGION_ID)

        BW         = imbinarize(double(Mice_STRUCT_T(mouseidx).ATLAS2COMET ==regionidx));
        boundaries = bwboundaries(BW);

        if (~isempty(boundaries))
            for k  = 1:1

                b  = boundaries{k};
                b(:,2) = b(:,2) + 985 - Bregma_Position(mouseidx,1);

                plot(b(:,2), b(:,1), 'Color' , [0 , 0 , 0]/255, 'LineWidth', 1); % Plot boundary with green line
                %     fill(b(:,2)+1, b(:,1)+1, Mice_STRUCT_T(1).COLOR(11,:),'FaceAlpha',0.5,'LineStyle',':');
                hold on

            end
            hold on
        end
        hold on
        regionidx
    end


    hold off

    %pause


end






close all

for learnidx=1:6


    figure(7000+learnidx)
    x              = Position_T(:,1);
    y              = Position_T(:,2);
    for regionidx=1:length(Mice_STRUCT_T(mouseidx).REGION_ID)

        BW         = imbinarize(double(Mice_STRUCT_T(mouseidx).ATLAS2COMET ==regionidx));
        boundaries = bwboundaries(BW);

        if (~isempty(boundaries))
            for k  = 1:1

                b  = boundaries{k};
                b(:,2) = b(:,2) + 985 - Bregma_Position(mouseidx,1);
                plot(b(:,2), b(:,1), 'Color' , [0 0 0]/255, 'LineWidth', 1); % Plot boundary with green line
                %     fill(b(:,2)+1, b(:,1)+1, Mice_STRUCT_T(1).COLOR(11,:),'FaceAlpha',0.5,'LineStyle',':');
                hold on

            end
            hold on
        end
        hold on
        regionidx
    end


    hold on


    nodeDegrees       = degree( G_neurons_learn{learnidx} );
    nodeDegrees_scale = rescale(nodeDegrees,0.1 , 1);
    scatter(x, y, 200*(nodeDegrees_scale), color_inuse_GROUP(learnidx,:), 'filled',  'MarkerFaceAlpha', 'flat', 'AlphaData', nodeDegrees_scale);

    axis equal
    %pause


end




Region_index_chosen = [1 1*regionidx]

    k            = 1;
    k2           = 1;
    vectors      = [];
    vectors2     = [];

    for mouseidx = 1


        Trace_Region9        = []
        Trace_Region10       = []
        Trace_Region9_learn  = []
        Trace_Region10_learn = []
        degree_centrality_9  = []
        degree_centrality_10 = []
        degree_centrality_T  = []

        regionidx1           = Region_index_chosen(1)%3%11%5%13%15%17%17%15


        Trace_Region9(:,:)   = Mice_STRUCT_T(mouseidx).NEURON.Trace_unify(Mice_STRUCT_T(mouseidx).Neuron_REGION_ID==CortexID_Mouse_Sorted(mouseidx,(regionidx1)) , logical(Frame_Inuse_Mouse{mouseidx}));
        Positions_9          = Mice_STRUCT_T(mouseidx).COORDINATES(Mice_STRUCT_T(mouseidx).Neuron_REGION_ID==CortexID_Mouse_Sorted(mouseidx,regionidx1),:);


        regionidx2           = Region_index_chosen(2)%4%12%6%14%16%18%18%16
        Trace_Region10(:,:)  = Mice_STRUCT_T(mouseidx).NEURON.Trace_unify(Mice_STRUCT_T(mouseidx).Neuron_REGION_ID==CortexID_Mouse_Sorted(mouseidx,regionidx2) , logical(Frame_Inuse_Mouse{mouseidx}));
        Positions_10         = Mice_STRUCT_T(mouseidx).COORDINATES(Mice_STRUCT_T(mouseidx).Neuron_REGION_ID==CortexID_Mouse_Sorted(mouseidx,regionidx2),:);


        Trace_Region_T       = [Trace_Region9;Trace_Region10];
        Position_T           = [Positions_9;Positions_10];


        x                    = Position_T(:,1);
        y                    = Position_T(:,2);

        for learnidx         = 2:max(Frame_chosen_index)

            Trace_Region_T_learn    = Trace_Region_T(:,Frame_chosen_index==learnidx);

            correlation_matrix_T    = corr(Trace_Region_T_learn');
            correlation_matrix_T(isnan(correlation_matrix_T)) = 0;

            binary_matrix_T         = correlation_matrix_T ;

            binary_matrix_T_adjcent = binary_matrix_T - diag(diag(binary_matrix_T));

            degree_centrality_T(:,learnidx) = (sum(binary_matrix_T_adjcent, 2))/(size(binary_matrix_T_adjcent,1)-1);;

            degree_centrality_in(regionidxL,learnidx)  = mean(sum(binary_matrix_T_adjcent(1:size(Positions_9,1),1:size(Positions_9,1)), 2))/size(Positions_9,1);
            degree_centrality_out(regionidxL,regionidx,learnidx) = mean(sum(binary_matrix_T_adjcent(1:size(Positions_9,1),1+size(Positions_9,1):end), 2))/size(Positions_9,1);

            binary_matrix_T(binary_matrix_T<0.8) = 0;

            [source, target, weight] = find(triu( binary_matrix_T ));
            G_neurons                = graph( source, target, weight );

            temp  = table2array(G_neurons.Edges);

%             for i = 1:size(temp, 1)
% 
% 
%                 if(temp(i,1) <= size(Positions_9,1)&temp(i,2)>size(Positions_9,1))
%                     points1   = [x(temp(i,1)) , y(temp(i,1))];
%                     points2   = [x(temp(i,2)) , y(temp(i,2))];
%                     vectors(k, :) = points2 - points1; % 边的向量
%                     k=k+1
%                 elseif(temp(i,1)>size(Positions_9,1)&temp(i,2)<=size(Positions_9,1))
%                     points2 = [x(temp(i,1)) , y(temp(i,1))];
%                     points1 = [x(temp(i,2)) , y(temp(i,2))];
%                     vectors(k, :) = points2 - points1; % 边的向量
%                     k=k+1
%                 end
% 
% 
%             end
% 
% 
% 
%             for i=1:size(temp, 1)
% 
%                 randomIndex     = randi(size(Position_T,1)); % 生成一个随机索引
%                 points1_rand(i) = randomIndex; % 选择对应的元素
%                 randomIndex     = randi(size(Position_T,1)); % 生成一个随机索引
%                 points2_rand(i) = randomIndex; % 选择对应的元素
% 
%             end
% 
%             for i = 1:size(temp, 1)
% 
% 
%                 if( points1_rand(i)<=size(Positions_9,1)&points2_rand(i)>size(Positions_9,1))
%                     points1 = [x(points1_rand(i)) , y(points1_rand(i))];
%                     points2 = [x(points2_rand(i)) , y(points2_rand(i))];
%                     vectors2(k2, :) = points2 - points1; % 边的向量
%                     k2=k2+1
%                 elseif( points1_rand(i)>size(Positions_9,1)&points2_rand(i)<=size(Positions_9,1))
%                     points2 = [x(points1_rand(i)) , y(points1_rand(i))];
%                     points1 = [x(points2_rand(i)) , y(points2_rand(i))];
%                     vectors2(k2, :) = points2 - points1; % 边的向量
%                     k2=k2+1
%                 end
% 
% 
%             end


            %         hold on
            % Region_index_chosen = SSP_bfd_index
            %     Trace_Region9        = []
            %     Trace_Region10       = []
            %     Trace_Region9_learn  = []
            %     Trace_Region10_learn = []
            %     degree_centrality_9  = []
            %     degree_centrality_10 = []
            %     degree_centrality_T  = []
            %
            %     regionidx = Region_index_chosen(1)%3%11%5%13%15%17%17%15
            %
            %     Trace_Region9(:,:)  = Mice_STRUCT_T(mouseidx).NEURON.Trace_unify(Mice_STRUCT_T(mouseidx).Neuron_REGION_ID==CortexID_Mouse_Sorted(mouseidx,regionidx) , logical(Frame_Inuse_Mouse{mouseidx}));
            %
            %     regionidx = Region_index_chosen(2)%4%12%6%14%16%18%18%16
            %
            %     Trace_Region10(:,:)  = Mice_STRUCT_T(mouseidx).NEURON.Trace_unify(Mice_STRUCT_T(mouseidx).Neuron_REGION_ID==CortexID_Mouse_Sorted(mouseidx,regionidx) , logical(Frame_Inuse_Mouse{mouseidx}));
            %
            %     Trace_Region_T = [Trace_Region9;Trace_Region10];
            %
            %   Trace_Region_T_learn  = Trace_Region_T(:,Frame_chosen_index==learnidx);
            %
            %         correlation_matrix_T  = corr(Trace_Region_T_learn');
            %         correlation_matrix_T(isnan(correlation_matrix_T))=0;
            %
            %         binary_matrix_T                   = correlation_matrix_T ;
            % binary_matrix_T(binary_matrix_T<0.6) = 0;
            %
            % [source, target, weight] = find(triu( binary_matrix_T ));
            % G_neurons = graph( source, target, weight );
            % p2=plot(G_neurons);%,'EdgeLabel', round(G.Edges.Weight,1)
            % % %pause

%             figure(500)
% 
%             p=plot(G_neurons, 'XData', x+ 985 - Bregma_Position(mouseidx,1), 'YData', y, 'EdgeAlpha', 0.2);%,'EdgeLabel', round(G.Edges.Weight,1)
%             %     axis equal
% %             title([Cortex_13regions_new_cell{regionidx} '-M: ' num2str(mouseidx) '-Learn: ' num2str(learnidx)])
%             hold on
% %             angles = atan2(vectors(:, 2), vectors(:, 1));
% %             angles_deg = abs(rad2deg(angles)); % 弧度转角度
% %             figure(888+0)
% %             histogram(angles_deg,'BinWidth', 2)
% % %             title(Cortex_13regions_new_cell{regionidx})
% % 
% %             angles = atan2(vectors2(:, 2), vectors2(:, 1));
% %             angles_deg = abs(rad2deg(angles)); % 弧度转角度
% %             figure(888+1)
% %             histogram(angles_deg,'BinWidth', 2)
% % %             title(Cortex_13regions_new_cell{regionidx})

%             %pause

        end



        %  %pause


        % binary_matrix_T_Mice(:,:,mouseidx) = binary_matrix_T;


        %     figure(500);
        %     hold on;
        %     for learnidx=1:max(Frame_chosen_index)
        %         plot(degree_centrality_9(:, learnidx), 'DisplayName', ['Stage ' num2str(learnidx)]);
        %     end
        %     xlabel('Neuron Index');
        %     ylabel('Degree Centrality');
        %     title('Degree Centrality Across Stages');
        %     legend show;
        %     hold off;
        %
        %     % 计算每个阶段的平均度中心性
        %     mean_degree_centrality_9(mouseidx,:) = sum(degree_centrality_9, 1);
        %     disp('Mean Degree Centrality for Each Stage:');
        %     disp(mean_degree_centrality_9);
        %
        %
        %     figure(501);
        %     hold on;
        %     for learnidx=1:max(Frame_chosen_index)
        %         plot(degree_centrality_10(:, learnidx), 'DisplayName', ['Stage ' num2str(learnidx)]);
        %     end
        %     xlabel('Neuron Index');
        %     ylabel('Degree Centrality');
        %     title('Degree Centrality Across Stages');
        %     legend show;
        %     hold off;
        %
        %     % 计算每个阶段的平均度中心性
        %     mean_degree_centrality_10(mouseidx,:) = sum(degree_centrality_10, 1);
        %     disp('Mean Degree Centrality for Each Stage:');
        %     disp(mean_degree_centrality_10);
        %
        %
        %     figure(502);
        %     hold on;
        %     for learnidx=1:max(Frame_chosen_index)
        %         plot(degree_centrality_T(:, learnidx), 'DisplayName', ['Stage ' num2str(learnidx)]);
        %     end
        %     xlabel('Neuron Index');
        %     ylabel('Degree Centrality');
        %     title('Degree Centrality Across Stages');
        %     legend show;
        %     hold off;

        % 计算每个阶段的平均度中心性
        %     mean_degree_centrality_T(mouseidx,:) = mean(degree_centrality_T, 1);
        %     disp('Mean Degree Centrality for Each Stage:');
        %     disp(degree_centrality_T);
        %     %pause

    end


% end
% 
% end




for learnidx=1:6

    figure(888+learnidx)
    imagesc(degree_centrality_out(:,:,learnidx))

end






for mouseidx   = 1:5


    degree_centrality_T  = []
    Trace_all            = []
    Trace_all            = Mice_STRUCT_T(mouseidx).NEURON.Trace_unify(: , logical(Frame_Inuse_Mouse{mouseidx}));

    Trace_Region_T_learn = []
    degree_centrality_T  = []
    correlation_matrix_T = []

    for learnidx=1:max(Frame_chosen_index)

        Trace_Region_T_learn = Trace_all(:,Frame_chosen_index==learnidx);

        correlation_matrix_T = corr(Trace_Region_T_learn');
        correlation_matrix_T(isnan(correlation_matrix_T))=0;

        binary_matrix_T = correlation_matrix_T ;

        degree_centrality_T(:, learnidx)   = (sum(binary_matrix_T, 2)-1)/(size(binary_matrix_T,1)-1);



%         n_random               = 10; % 随机网络的数量
%         binary_matrix_T_adjcent = binary_matrix_T - diag(diag(binary_matrix_T));
%         binary_matrix_T_adjcent = binary_matrix_T_adjcent>0.5;
%         [gamma, lambda, sigma] = small_world_properties(binary_matrix_T_adjcent ,n_random);
% 
% 
%         gamma_mice_cortex( mouseidx,learnidx)  =    gamma;
%         lambda_mice_cortex(mouseidx,learnidx) =    lambda;
%         sigma_mice_cortex( mouseidx,learnidx)  =    sigma;



    end


    mean_degree_centrality_T(mouseidx,:) = mean(degree_centrality_T, 1);
    mouseidx


end







aa=mean_degree_centrality_T'

degree_region_centrality_T=[]
degree_region_centrality_T_mice=[]
correlation_matrix_T_Mice_Learn=[]
close all

for mouseidx   = 1:5


    Trace_unify_Region = Mice_STRUCT_T(mouseidx).NEURON.Trace_unify_Region_Avg(CortexID_Mouse_Sorted(mouseidx,:) , logical(Frame_Inuse_Mouse{mouseidx}));



    for learnidx=1:6


        Trace_Region_T_learn = Trace_unify_Region(:,Frame_Learn_index==learnidx);

        correlation_matrix_T = corr(Trace_Region_T_learn');

        threshold = 0.8;

        binary_matrix_T = correlation_matrix_T% > threshold;

        degree_region_centrality_T(:, learnidx) = (sum(binary_matrix_T, 2)-1)/(size(binary_matrix_T,1)-1);

        degree_region_centrality_T_mice(:,mouseidx,learnidx) = degree_region_centrality_T(:, learnidx);
      
        figure(500+learnidx)
        imagesc(binary_matrix_T)

         


    end
%     %pause
    sum_degree_regions_centrality_T(mouseidx,:) = sum(degree_region_centrality_T, 1);

end



Cortex_Incommon_New = Cortex_Incommon([3:4, 1:2, 11:12, 15:16, 13:14, 19:20,    17:18,9:10, 7:8, 5:6,  21:22, 25:26, 23:24]);

Cortex_Incommon_New = Cortex_Incommon([3:4, 1:2, 11:12, 15:16, 13:14, 19:20,17:18,9:10,     7:8, 5:6,  21:22, 25:26, 23:24]);

Cortex_Incommon_New = Cortex_Incommon([3:4, 1:2, 15:16,  11:12,  13:14, 17:18,9:10, 19:20,7:8,   5:6,  21:22, 25:26, 23:24]);

Cortex_Incommon_New = Cortex_Incommon([3:4, 1:2, 15:16,  11:12,  13:14, 17:18,9:10, 19:20,    21:22, 25:26, 23:24,5:6, 7:8]);

% Cortex_Incommon_New = Cortex_Incommon([3:4, 1:2, 15:16,  11:12,  13:14, 17:18,9:10, 19:20, 27,   21:22, 25:26, 23:24,5:6, 7:8]);


[~, indexInArray2] = ismember(Cortex_Incommon_New,Cortex_Incommon_Sorted);

% 输出结果
disp("strArray1 中元素在 strArray2 中的索引：");
disp(indexInArray2);



Frame_chosen_index              = Frame_Learn_index
sum_degree_regions_centrality_T = []
correlation_matrix_T_Mice_Learn = []
degree_region_centrality_T_mice = []
binary_matrix_T                 = []

for mouseidx    = 1:5

    Trace_Region_Avg_Tmouse = Trace_Region_Avg(:,:,mouseidx) ;

     Trace_Region_Avg_Tmouse_group = Trace_Region_Avg_Tmouse(indexInArray2, :);

    Trace_Region_Avg_TmouseAvg_chosen=[]

    for regionidx = 1:8

        Trace_Region_Avg_TmouseAvg_chosen(regionidx,:) = mean(Trace_Region_Avg_Tmouse_group([2*regionidx-1 2*regionidx],:),1);

    end
% regionidx=9;
%      Trace_Region_Avg_TmouseAvg_chosen(regionidx,:) = mean(Trace_Region_Avg_Tmouse_group([2*regionidx-1],:),1);

   for regionidx = 9:13

        Trace_Region_Avg_TmouseAvg_chosen(regionidx,:) = mean(Trace_Region_Avg_Tmouse_group([2*regionidx-1 2*regionidx],:),1);

    end

    degree_region_centrality_T  = []

    for learnidx=1:6


        Trace_Region_T_learn = Trace_Region_Avg_TmouseAvg_chosen(:,Frame_chosen_index==learnidx);

        correlation_matrix_T = corr(Trace_Region_T_learn');

        threshold = 0;

        %         binary_matrix_T = correlation_matrix_T > threshold;
        binary_matrix_T(:,:,mouseidx) = correlation_matrix_T ;

        degree_region_centrality_T(:, learnidx) = (sum( binary_matrix_T(:,:,mouseidx) , 2)-1)/(size( binary_matrix_T(:,:,mouseidx) ,1)-1);

        degree_region_centrality_T_mice(:,mouseidx,learnidx) = degree_region_centrality_T(:, learnidx);

        figure(500+learnidx)
        imagesc( binary_matrix_T(:,:,mouseidx) )

        correlation_matrix_T_Mice_Learn(:,:,learnidx,mouseidx) = correlation_matrix_T;

    end


    %pause

    sum_degree_regions_centrality_T(mouseidx,:) = sum(degree_region_centrality_T, 1);


end






close all

for learnidx=1:6
        figure(7000+learnidx)
        temp = mean(correlation_matrix_T_Mice_Learn(:,:,learnidx,:),4);
        colormap(parula)
%  colormap(jet)
        imagesc( temp >0.52)%>0.80
        axis square;
        xticks(1:13); % 设置 X 轴刻度
yticks(1:13); % 设置 Y 轴刻度
xticklabels(  [Cortex_Incommon_New(1:2:17); Cortex_Incommon_New(18:2:26) ]  ); % 自定义 X 轴刻度标签
yticklabels(  [Cortex_Incommon_New(1:2:17); Cortex_Incommon_New(18:2:26) ]  ); % 自定义 Y 轴刻度标签

end



correlation_matrix_T_Mice_Learn_avg = mean(  correlation_matrix_T_Mice_Learn , 4);
degree_MOs2Comparator=[]
degree_MOs2Integrator=[]
degree_Comparator2Integrator=[]
degree_Comparator2VISa=[]
degree_Comparator2RSPd=[]
degree_Comparator2SSpbfd=[]
for learnidx=1:6

    degree_MOs2Comparator(learnidx,:) = mean(  correlation_matrix_T_Mice_Learn(1, 2:8,learnidx,:) ,2 ) /mean(correlation_matrix_T_Mice_Learn(:,:,learnidx,:),'all')
    degree_MOs2Integrator(learnidx,:) = mean(  correlation_matrix_T_Mice_Learn(1, 9:11,learnidx,:)  ,2 ) /mean(correlation_matrix_T_Mice_Learn(:,:,learnidx,:),'all')
    degree_Comparator2Integrator(learnidx,:) = mean( mean( correlation_matrix_T_Mice_Learn(2:8, 9:11,learnidx,:) ,2  ),1  ) /mean(correlation_matrix_T_Mice_Learn(:,:,learnidx,:),'all')
    degree_Comparator2VISa(learnidx,:) = mean( mean( correlation_matrix_T_Mice_Learn(2:7, 8,learnidx,:) ,2  ),1  )/mean(correlation_matrix_T_Mice_Learn(:,:,learnidx,:),'all')
    degree_Comparator2RSPd(learnidx,:) = mean( mean( correlation_matrix_T_Mice_Learn(2:8, 13,learnidx,:)  ,2  ),1  )/mean(correlation_matrix_T_Mice_Learn(:,:,learnidx,:),'all')
    degree_Comparator2SSpbfd(learnidx,:) = mean( mean( correlation_matrix_T_Mice_Learn(2:8, 12,learnidx,:) ,2  ),1  )/mean(correlation_matrix_T_Mice_Learn(:,:,learnidx,:),'all')

    degree_VISa2MOs(learnidx,:) = mean(  correlation_matrix_T_Mice_Learn(8, 1,learnidx,:) ,2 ) /1
    degree_VISa2Integrator(learnidx,:) = mean(  correlation_matrix_T_Mice_Learn(8, 9:11,learnidx,:)  ,2 ) /1
    degree_VISa2Comparator(learnidx,:) = mean( mean( correlation_matrix_T_Mice_Learn(8, 2:7,learnidx,:) ,2  ),1  ) /1
  
    degree_VISa2RSPd(learnidx,:) = mean( mean( correlation_matrix_T_Mice_Learn(8, 13,learnidx,:)  ,2  ),1  )/1
    degree_VISa2SSpbfd(learnidx,:) = mean( mean( correlation_matrix_T_Mice_Learn(8, 12,learnidx,:) ,2  ),1  )/1

end


nodeIds = [1; 2; 3; 4;5;6;7;8;9;10;11;12;13]; % 节点 ID
nodeLabels = Cortex_13regions_new_cell; % 节点标签
nodeCategories = {'Group1'; 'Group2'; 'Group2'; 'Group2'; 'Group2'; 'Group2'; 'Group2'; 'Group2';'Group3';'Group3';'Group3';'Group4';'Group5'}; % 节点类别

color_inuse_GROUP                  =[]
color_inuse_GROUP                  = color_inuse([19 3 2 16 24],:);
color_inuse_GROUP(6,:)             = [0.5 0.5 0.5];



correlation_matrix_T_Mice_Learn_avg = mean(correlation_matrix_T_Mice_Learn,4);


close all

for learnidx=1:6


    temp = correlation_matrix_T_Mice_Learn_avg(:,:,learnidx)  - diag(diag(correlation_matrix_T_Mice_Learn_avg(:,:,learnidx) ));
  
 temp = correlation_matrix_T_Mice_Learn_avg(:,:,learnidx)  ;


    Degree_13regions_T(:, learnidx) = (sum(temp , 2));
%     temp(temp<0.8) = 0;
temp(temp<0.8) = 0;
    correlation_matrix_T_Mice_Learn_avg_GRAPH(:,:,learnidx)  = (temp);
    % 获取边的信息
    [source, target, weight] = find(triu( correlation_matrix_T_Mice_Learn_avg_GRAPH(:,:,learnidx) ));

    % 创建边列表
    edgeList = [source, target, weight];

    folder2 = 'E:\桌面工作\1-郭长亮项目文件夹\1-项目1-CM2SCOPE\1-V1-CM2scope\5-论文撰写\Manuscript\Figures\Figure4\';
    % 保存为 CSV 文件
    csvwrite([folder2 'learn-' num2str(learnidx) '-edge_list.csv'], edgeList);

    % 添加表头
    fid = fopen([folder2 'learn-' num2str(learnidx) '-edge_list.csv'], 'w');
    fprintf(fid, 'Source,Target,Weight\n');
    fclose(fid);

    % 追加数据
    dlmwrite([folder2 'learn-' num2str(learnidx) '-edge_list.csv'], edgeList, '-append', 'delimiter', ',');

    G = graph( source, target, weight, Cortex_13regions_new_cell );

    % y = [ 5  5  4  4  4  3  3  2  1  1  1  2  2]; % 节点的 x 坐标
    % x = [ 2  4  1  3  5  2  4  1  1  3  5  3  5]; % 节点的 y 坐标

    y = [ 5  5  4  3.5  4  2.6  2.6  1.6  1.2  1  1  2  2]; % 节点的 x 坐标
    x = [ 2  4  1  3    5  2    4    3    3    2  4  1  5]; % 节点的 y 坐标

    y = [ 6  5  4  3.5  4  2    2.6  1.6  1.2  1  1  2  2]; % 节点的 x 坐标
    x = [ 3  3  1  3    5  3    3    3    3    2  4  1  5]; % 节点的 y 坐标

    y = [ 6  5    3.5  4  4  2.5  3.5  1.5   1.2  1  1  2  2]; % 节点的 x 坐标
    x = [ 3  2.8  2    1  5  3.3  4    2.8   3.2  2  4  1  5]; % 节点的 y 坐标

    
    nodeDegrees = degree(G);
    nodeDegrees_region_learn(:,learnidx) = nodeDegrees;
    nodesize = floor(rescale(nodeDegrees)*10+1);
    % 可视化图
    handle2 = figure(900+learnidx);
    handle2.Position = [100 100 500 600];
    p=plot(G, 'Layout', 'force', 'NodeLabel',  G.Nodes.Name);%,'EdgeLabel', round(G.Edges.Weight,1)
    layout(p,'force','UseGravity',true)
    p.NodeFontSize=12
    p.XData=x;
    p.YData=y;
    p.NodeColor=[color_inuse_GROUP(1,:);
        color_inuse_GROUP(4,:);
        color_inuse_GROUP(4,:);
        color_inuse_GROUP(4,:);
        color_inuse_GROUP(4,:);
        color_inuse_GROUP(4,:);
        color_inuse_GROUP(4,:);
        color_inuse_GROUP(4,:);
        color_inuse_GROUP(3,:);
        color_inuse_GROUP(3,:);
        color_inuse_GROUP(3,:);
        color_inuse_GROUP(2,:);
        color_inuse_GROUP(5,:)
        ];
    p.MarkerSize=nodesize;
    p.EdgeColor='k';
    p.LineWidth=weight*2;
    % plot(G, 'NodeFontSize', 12,'XData', x, 'YData', y, 'NodeColor', 'k', 'MarkerSize', nodesize, 'EdgeColor', 'k', 'LineWidth', 1,'NodeLabel',  G.Nodes.Name);
    % title('包含孤立节点的图');
    set(gca, 'XColor', 'none');
    set(gca, 'XTickLabel', []);
    set(gca, 'YColor', 'none');
    set(gca, 'YTickLabel', []);




handle3=figure(1000+learnidx);
handle3.Position = [100 100 500 600];
hold on;

nodes = [x;y]'
% 绘制曲线边
for i = 1:size(source, 1)
    % 获取边的起点和终点
    startNode = nodes(source(i), :);
    endNode   = nodes(target(i), :);
    
    % 计算控制点（用于贝塞尔曲线）
    controlPoint = (startNode + endNode) / 2 + [0, 0.5];  % 调整控制点位置
    
    % 生成贝塞尔曲线
    t = linspace(0, 1, 100);
    t = t';

    curve = (1 - t).^2 .* repmat(startNode,length(t),1) + 2 * (1 - t) .* t .* repmat(controlPoint,length(t),1) + t.^2 .* repmat(endNode,length(t),1);
        % 绘制曲线
        
    plot(curve(:, 1), curve(:, 2), 'k', 'LineWidth', weight(i)*2);

end
hold on;
scatter(x, y, nodesize*15, [color_inuse_GROUP(1,:);
        color_inuse_GROUP(4,:);
        color_inuse_GROUP(4,:);
        color_inuse_GROUP(4,:);
        color_inuse_GROUP(4,:);
        color_inuse_GROUP(4,:);
        color_inuse_GROUP(4,:);
        color_inuse_GROUP(4,:);
        color_inuse_GROUP(3,:);
        color_inuse_GROUP(3,:);
        color_inuse_GROUP(3,:);
        color_inuse_GROUP(2,:);
        color_inuse_GROUP(5,:)
        ] , 'filled');
% % 设置图形属性
% title('手动绘制贝塞尔曲线边');
% axis equal;
% hold off;
 xlim([0 6])
  ylim([0 7])
   set(gca, 'XColor', 'none');
    set(gca, 'XTickLabel', []);
    set(gca, 'YColor', 'none');
    set(gca, 'YTickLabel', []);

end



% 创建图


nodeIds = [1; 2; 3; 4;5;6;7;8;9;10;11;12;13]; % 节点 ID
nodeLabels = Cortex_13regions_new_cell; % 节点标签
nodeCategories = {'Group1'; 'Group2'; 'Group2'; 'Group2'; 'Group2'; 'Group2'; 'Group2'; 'Group2';'Group3';'Group3';'Group3';'Group4';'Group5'}; % 节点类别

% 创建节点表
nodeTable = table(nodeIds, nodeLabels, nodeCategories, ...
    'VariableNames', {'Id', 'Label', 'Category'});

% 保存节点表为 CSV 文件
writetable(nodeTable, [folder2 'nodes_list.csv']);



correlation_matrix_T_Mice_Learn_avg





%%========================Network analysis=================================



Correlation_regions_mice=[]

for mouseidx=1:5

    Trace_Region_Avg_Tmouse = Trace_Region_Avg(:,:,mouseidx) ;

    for regionidx = 1:13


        Trace_Region_Avg_TmouseAvg_chosen(regionidx,:) = mean( Trace_Region_Avg_Tmouse([2*regionidx-1 2*regionidx],:) , 1);
        Trace_Region_Avg_TmouseAvg_avg(regionidx,:) = mean( Trace_Region_Avg_Tmouse_avg([2*regionidx-1 2*regionidx],:) , 1);

    end


    for learnidx=1:6


        temp =  corr(   (  Trace_Region_Avg_TmouseAvg_chosen(:,Frame_Learn_index==learnidx)  )'   );
        Correlation_regions_mice(:,:,learnidx,mouseidx) = temp;

        temp =  corr(   (  Trace_Region_Avg_TmouseAvg_avg(:,Frame_Learn_index==learnidx)  )'   );
        Correlation_regions_mice_avg(:,:,learnidx) = temp;


        temp =  pdist2(   (  Trace_Region_Avg_TmouseAvg_chosen(:,Frame_Learn_index==learnidx)  )  , (  Trace_Region_Avg_TmouseAvg_chosen(:,Frame_Learn_index==learnidx)  ) );
        pdist_regions_mice(:,:,learnidx,mouseidx) = temp;

        temp =  pdist2(   (  Trace_Region_Avg_TmouseAvg_avg(:,Frame_Learn_index==learnidx)  ) , (  Trace_Region_Avg_TmouseAvg_chosen(:,Frame_Learn_index==learnidx)  )  );
        pdist_regions_mice_avg(:,:,learnidx) = temp;

% 
%            temp =  corr(   (  Trace_Region_Avg_TmouseAvg_chosen(:,Frame_Learn_index==learnidx)  )  , (  Trace_Region_Avg_TmouseAvg_chosen(:,Frame_Learn_index==learnidx)  ) , 'Type', 'Kendall');
%         Kendall_regions_mice(:,:,learnidx,mouseidx) = temp;
% 
%         temp =  pdist2(   (  Trace_Region_Avg_TmouseAvg_avg(:,Frame_Learn_index==learnidx)  ) , (  Trace_Region_Avg_TmouseAvg_chosen(:,Frame_Learn_index==learnidx)  ) , 'Type', 'Kendall' );
%         Kendall_regions_mice_avg(:,:,learnidx) = temp;


    end


end



region_reindex = [6 2 4 11 12 10 13 9 3 7 1 5 8]

for mouseidx=1:5
    Correlation_regions_mice_chosen = Correlation_regions_mice(:,:,:,mouseidx);

    for learnidx=1:6

        figure(1000+learnidx)
        temp = Correlation_regions_mice_chosen(:,:,learnidx);
        temp_reindex = temp(region_reindex,region_reindex);
        imagesc(temp_reindex); % 绘制热图
        colorbar; % 添加颜色条

    end
    %pause
end

    for learnidx=1:6

        figure(2000+learnidx)
        temp = Correlation_regions_mice_avg(:,:,learnidx);
        temp_reindex = temp(region_reindex,region_reindex);
        imagesc(temp_reindex); % 绘制热图
        colorbar; % 添加颜色条

    end

        for learnidx=1:6

        figure(3000+learnidx)
        temp = pdist_regions_mice_avg(:,:,learnidx);
        temp_reindex = temp(region_reindex,region_reindex);
        imagesc(temp_reindex); % 绘制热图
        colorbar; % 添加颜色条

    end


Correlation_regions_mice_avg_thr = (Correlation_regions_mice_avg(:,:,6)-Correlation_regions_mice_avg(:,:,1))
figure(1002)
imagesc( Correlation_regions_mice_avg_thr ); % 绘制热图
colorbar; % 添加颜色条



 for regionidx = 1:13


     Trace_Region_Avg_TmouseAvg_chosen = Trace_Region_Avg_TmouseAvg([2*regionidx-1 2*regionidx],:);
     x  = 1:1:length(Trace_Region_Avg_TmouseAvg); % x 轴
     y1 = mean( Trace_Region_Avg_TmouseAvg , 1); % 曲线 1
     y2 = mean(Trace_Region_Avg_TmouseAvg_chosen,1); % 曲线 2（带噪声）

     % 参数设置
window_size = 400; % 窗口大小
step = 1; % 滑动步长
num_windows = floor((length(x) - window_size) / step) + 1; % 窗口数量
     % 初始化存储差异的数组
     mse_values = zeros(1, num_windows); % 存储每个窗口的 MSE
  

     % 滑动窗口计算差异
     for i = 1:num_windows
         % 获取当前窗口的数据
         start_idx = (i-1) * step + 1;
         end_idx = start_idx + window_size - 1;
         window_y1 = y1(start_idx:end_idx);
         window_y2 = y2(start_idx:end_idx);

         % 计算当前窗口的 MSE 和 MAE
         mse_values(i) = mean((window_y2 - window_y1)); % MSE
     
    
    end

     % 可视化结果
     window_centers = x(1:step:end-window_size+1) + (window_size/2) * (x(2)-x(1)); % 窗口中心位置

     figure( 900 + regionidx );
     subplot(2, 1, 1);
     plot(x, y1, 'b-', 'LineWidth', 2); % 绘制曲线 1
     hold on;
     plot(x, y2, 'r--', 'LineWidth', 2); % 绘制曲线 2
     title('原始曲线');
     legend('曲线 1', '曲线 2');
     xlabel('X');
     ylabel('Y');

     subplot(2, 1, 2);
     plot(window_centers, mse_values, 'g-', 'LineWidth', 2); % 绘制 MSE
   
     title('滑动窗口差异');
     legend('MSE', 'MAE');
     xlabel('X');
     ylabel('差异值');

    t1 = title(Regions_overlap_inseq{2*regionidx-1})
 end







 regionidx        = 27

 handle           = figure(701+regionidx)
 handle.Position  = [100 100 800 200];

 CueShock_marker_color_Y(miniscopecuebeginFrame, miniscopecueendFrame,miniscopeshockbeginFrame,miniscopeshockendFrame,0.15);

 plot(  mean( Trace_Region_Avg_TmouseAvg , 1)   , 'LineWidth',1,'Color', 'k'  )
 hold on

 Trace_Region_Avg_TmouseAvg_chosen = Trace_Region_Avg_TmouseAvg(regionidx,:);
 plot(  mean(Trace_Region_Avg_TmouseAvg_chosen,1)   , 'LineWidth',1,'Color', 'r'  )

 t1              = title(Regions_overlap_inseq{regionidx})


 Trace_unify_Region_Avg_VISrlLeft  = Mice_STRUCT_T(mouseidx).NEURON.Trace_unify_Region_Avg(82,logical(Frame_Inuse_Mouse{mouseidx}));  

 handle          = figure(750)
 handle.Position = [100 100 800 200];

 CueShock_marker_color_Y(miniscopecuebeginFrame, miniscopecueendFrame,miniscopeshockbeginFrame,miniscopeshockendFrame,0.15);

 plot(  mean( Trace_Region_Avg_TmouseAvg , 1)   , 'LineWidth',1,'Color', 'k'  )
 hold on
 plot(  mean(Trace_unify_Region_Avg_VISrlLeft,1)   , 'LineWidth',1,'Color', 'r'  )







for regionidx=1:27

    writematrix(squeeze( Trace_Region_Avg(regionidx,:,:) ) , ['E:\桌面工作\1-郭长亮项目文件夹\1-项目1-CM2SCOPE\1-V1-CM2scope\5-论文撰写\Manuscript\Figures\Figure3\' 'Trace_Region_Avg.xls'] ,'Sheet',regionidx);

    regionidx

end



Non_Responder_Regions              = [  11 12          ];
Average_Regions                    = [  13 14                 ];
Integrator_Regions                 = [  9 10 15 16  1 2              ];
Comparison_Regions                 = [  3  4 7  8 17 18 19 20 21 22  23 24 25 26 27 ];
RSPd_Regions                       = [  5  6                      ];


% Non_Responder_Regions              = [   3 4 11 12          ];
% Average_Regions                    = [ 13 14  25 26 27               ];
% Integrator_Regions                 = [  9 10 15 16  1 2              ];
% Comparison_Regions                 = [   7  8 17 18 19 20 21 22  23 24 ];
% RSPd_Regions                       = [  5  6                      ];
% 

Regions_GROUPS                     = {};
Regions_GROUPS{1}  = Non_Responder_Regions
Regions_GROUPS{2}  = Average_Regions
Regions_GROUPS{3}  = Integrator_Regions
Regions_GROUPS{4}  = Comparison_Regions
Regions_GROUPS{5}  = RSPd_Regions

Regions_GROUPS_NAMES = {};
Regions_GROUPS_NAMES{1}  = 'Non_Responder Regions'
Regions_GROUPS_NAMES{2}  = 'Average Regions'
Regions_GROUPS_NAMES{3}  = 'Integrator Regions'
Regions_GROUPS_NAMES{4}  = 'Comparison Regions'
Regions_GROUPS_NAMES{5}  = 'RSPd Regions'

% Other_Regions         = [  1  2                      ];
 Trace_Region_Avg_TmouseAvg = mean( Trace_Region_Avg , 3 );


 for regionidx = 1:5

     handle = figure(750+regionidx)
     handle.Position = [100 100 800 200];

     CueShock_marker_color_Y(miniscopecuebeginFrame, miniscopecueendFrame,miniscopeshockbeginFrame,miniscopeshockendFrame,0.15);

     plot(  mean( Trace_Region_Avg_TmouseAvg , 1)   , 'LineWidth',1,'Color', 'k'  )
     hold on

     Trace_Region_Avg_TmouseAvg_chosen = Trace_Region_Avg_TmouseAvg(Regions_GROUPS{regionidx} ,:);
     plot(  mean(Trace_Region_Avg_TmouseAvg_chosen,1)   , 'LineWidth',1,'Color', 'r'  )

     ylim([0.04 0.13])
     t1 = title(Regions_GROUPS_NAMES{regionidx})
set(gca, 'XColor', 'none');
     set(gca, 'XTickLabel', []);
     set(gca, 'YColor', 'none');
     set(gca, 'YTickLabel', []);

%      %pause

%      close(handle)

 end


[~, indices_Non_Responder_Regions] = ismember(  Regions_overlap_inseq(Non_Responder_Regions)  ,  Mice_STRUCT_T(1).REGION_NAME  );
[~, indices_Average_Regions]       = ismember(  Regions_overlap_inseq(Average_Regions)        ,  Mice_STRUCT_T(1).REGION_NAME  );
[~, indices_Integrator_Regions]    = ismember(  Regions_overlap_inseq(Integrator_Regions)     ,  Mice_STRUCT_T(1).REGION_NAME  );
[~, indices_Comparison_Regions]    = ismember(  Regions_overlap_inseq(Comparison_Regions)     ,  Mice_STRUCT_T(1).REGION_NAME  );
[~, indices_RSPd_Regions]          = ismember(  Regions_overlap_inseq(RSPd_Regions)           ,  Mice_STRUCT_T(1).REGION_NAME  );
% [~, indices_Other_Regions]         = ismember(  Regions_overlap_inseq(Other_Regions)          ,  Mice_STRUCT_T(1).REGION_NAME  );


color_inuse_GROUP                  =[]
color_inuse_GROUP                  = color_inuse([19 3 2 16 24],:);
color_inuse_GROUP(6,:)             = [0.5 0.5 0.5];

mouseidx   = 1

Neuron_REGION_ID_GROUP             = 6*ones(size(Mice_STRUCT_T(mouseidx).COORDINATES,1) , 1  );

for i = 1:length(indices_Non_Responder_Regions)

    Neuron_REGION_ID_GROUP(Mice_STRUCT_T(mouseidx).Neuron_REGION_ID == indices_Non_Responder_Regions(i)) = 1;

end

for i = 1:length(indices_Average_Regions)

    Neuron_REGION_ID_GROUP(Mice_STRUCT_T(mouseidx).Neuron_REGION_ID == indices_Average_Regions(i)) = 2;

end

for i = 1:length(indices_Integrator_Regions)

    Neuron_REGION_ID_GROUP(Mice_STRUCT_T(mouseidx).Neuron_REGION_ID == indices_Integrator_Regions(i)) = 3;

end

for i = 1:length(indices_Comparison_Regions)

    Neuron_REGION_ID_GROUP(Mice_STRUCT_T(mouseidx).Neuron_REGION_ID == indices_Comparison_Regions(i)) = 4;

end

for i = 1:length(indices_RSPd_Regions)

    Neuron_REGION_ID_GROUP(Mice_STRUCT_T(mouseidx).Neuron_REGION_ID == indices_RSPd_Regions(i)) = 5;

end


% 
% for i = 1:length(indices_Other_Regions)
% 
%     Neuron_REGION_ID_GROUP(Mice_STRUCT_T(mouseidx).Neuron_REGION_ID == indices_RSPd_Regions(i)) = 6;
% 
% end


figure(700)
    imshow(uint8(1* ones(size(Mice_STRUCT_T(mouseidx).ATLAS2COMET,1),size(Mice_STRUCT_T(mouseidx).ATLAS2COMET,2)) ));
    hold on;
    for regionidx=1:length(Mice_STRUCT_T(mouseidx).REGION_ID)

        BW = imbinarize(double(Mice_STRUCT_T(mouseidx).ATLAS2COMET ==regionidx));
        boundaries = bwboundaries(BW);

        if (~isempty(boundaries))
            for k = 1:1

                b = boundaries{k};
                plot(b(:,2), b(:,1), 'Color' , [0 1 1], 'LineWidth', 1); % Plot boundary with green line
                %     fill(b(:,2)+1, b(:,1)+1, Mice_STRUCT_T(1).COLOR(11,:),'FaceAlpha',0.5,'LineStyle',':');
                hold on

            end
            hold on
        end
        hold on
        regionidx
    end

    hold on
    scatter( Mice_STRUCT_T(mouseidx).COORDINATES(:,1) , Mice_STRUCT_T(mouseidx).COORDINATES(:,2),8, color_inuse_GROUP(Neuron_REGION_ID_GROUP,:) ,'filled' )
    hold off


   
    FIG_COLORBAR=figure(750);
    FIG_COLORBAR.Position = [500,50,180,900]
    COLOR_MOUSE = Mice_STRUCT_T(mouseidx).COLOR(unique(Mice_STRUCT_T(mouseidx).Neuron_REGION_ID) , :);
    scatter( repmat(1,1,size(COLOR_MOUSE,1)), 1:size(COLOR_MOUSE,1) , 100,COLOR_MOUSE ,'filled')

    text(repmat(1.15,1,size(COLOR_MOUSE,1)),[1:size(COLOR_MOUSE,1)]-0.1,Mice_STRUCT_T(mouseidx).REGION_NAME(unique(Mice_STRUCT_T(mouseidx).Neuron_REGION_ID)) , 'FontSize', 12 )
    axis off;


figure(799)
colormap(color_inuse); % Use the 'jet' colormap
colorbar;


%%=========================================================================
%%=========================================================================
%%=======================SHOCK=============================================
%%=======================SHOCK=============================================


close all

Activity_SHOCK_avg    = []
Activity_SHOCK_trial  = []
Trace_SHOCK_avg       = []
Zeropad_left          = 5;%10;
Zeropad_right         = 50;%50;

for mouseidx = 1:5

    for regionidx      = 1:27

        Activity_SHOCK = NaN( 20+Zeropad_left+Zeropad_right , 5 );
        Trace_SHOCK    = NaN(20+Zeropad_left+Zeropad_right,5);
        SHOCK_frame_z_trace = zeros(20+Zeropad_left+Zeropad_right,5);

        for trialidx=1:5


            SHOCK_frame   = find(Mice_STRUCT_T(mouseidx).NEURON.SHOCK(logical(Frame_Inuse_Mouse{mouseidx}))==trialidx);

            SHOCK_frame_z = [SHOCK_frame(1)-Zeropad_left:SHOCK_frame(1)-1, SHOCK_frame, SHOCK_frame(end)+1:SHOCK_frame(end)+Zeropad_right];

            temp          = Mice_STRUCT_T(mouseidx).NEURON.SHOCK(logical(Frame_Inuse_Mouse{mouseidx}));

            SHOCK_frame_z_trace(1:length(SHOCK_frame_z),trialidx) = .2*(temp(SHOCK_frame_z)>0);

            Activity_SHOCK(1:length(SHOCK_frame_z),trialidx) = single(Trace_Region_Avg(regionidx , SHOCK_frame_z , mouseidx))';

            Trace_SHOCK(1:length(SHOCK_frame_z),trialidx)    = single(Trace_unify_AVG(  mouseidx , SHOCK_frame_z  ))';
  
       
            %         pause

        end

            Activity_SHOCK_trial(:,:,regionidx,mouseidx) = Activity_SHOCK;

            Activity_SHOCK_avg(:,regionidx,mouseidx) = nanmean(Activity_SHOCK,2);
            Trace_SHOCK_avg(:,:,mouseidx) = Trace_SHOCK;

    end

% % %         figure(500)
% % % 
% % %         subplot(5,6,regionidx)
% % % 
% % %         title(Cortex_Incommon_Sorted{regionidx})
% % %         plot(Activity_SHOCK, 'LineWidth',2 ,'Color','r')
% % %         hold on
% % %         plot(Trace_SHOCK , 'LineWidth',2,'Color', 'k')
% % % 
% % %         title([Mice_STRUCT_T(mouseidx).NAME(4:end) '  ' Cortex_Incommon_Sorted{regionidx}])
% % % 
% % %         axis tight
% % %         %     axis off
% % %         ylim([0 0.2])
% % %         yticks([0:0.1:0.2])
% % %         %             ylim([-3 3])
% % %         %         yticks([-3:1:3])
% % % 
% % %         hold off

    end

%     pause
% end


Activity_SHOCK_trial_MOTION_M94 = reshape( Activity_SHOCK_trial(:,:,index_2_Regroup(1:4),1:5), [20+Zeropad_left+Zeropad_right, 5*4*5] );

Trace_SHOCK_avg_Whole_2D        = reshape( Trace_SHOCK_avg, [20+Zeropad_left+Zeropad_right, 5*5] );

Activity_SHOCK_trial_MOTION_SSP = reshape( Activity_SHOCK_trial(:,:,index_2_Regroup(5:14),1:5), [20+Zeropad_left+Zeropad_right, 5*10*5] );
Activity_SHOCK_trial_MOTION_RSP = reshape( Activity_SHOCK_trial(:,:,index_2_Regroup(15:16),1:5), [20+Zeropad_left+Zeropad_right, 5*2*5] );
Activity_SHOCK_trial_MOTION_VIS = reshape( Activity_SHOCK_trial(:,:,index_2_Regroup(21:26),1:5), [20+Zeropad_left+Zeropad_right, 5*6*5] );
Activity_SHOCK_trial_MOTION_VIS_w = reshape( Activity_SHOCK_trial(:,:,index_2_Regroup(19:27),1:5), [20+Zeropad_left+Zeropad_right, 5*9*5] );



Color_trial  = [255 51  51 ;
                255 153 51 ;
                255 102 255;
                0   102 0  ;
                0   0   255]/255


%------------------------------------------
for mouseidx=1:5

    figure( 600  + mouseidx )
    t1           = find(diff(SHOCK_frame_z_trace(:,1))>0)+1; % 开始时间
    t2           = find(diff(SHOCK_frame_z_trace(:,1))<0); % 结束时间
    % 添加阴影
    px           = [t1, t2, t2, t1]; % 阴影的 x 坐标
    py           = [0,  0,  1,  1 ]; % 阴影的 y 坐标
    patch( px, py, [178 102 255]/255, 'FaceAlpha', 1, 'EdgeColor', 'none' ); % 添加阴影
    hold on
    for trialidx = 1:5

        plot(  Trace_SHOCK_avg(:,trialidx,mouseidx)*7,'LineWidth', 7, 'Color',Color_trial(trialidx,:)  )
        % pause
        hold on

    end
    axis off;
end

%------------------------------------------
%------------------------------------------


close all
for mouseidx=1:5

    figure( 606  + mouseidx )
    t1           = find(diff(SHOCK_frame_z_trace(:,1))>0)+1; % 开始时间
    t2           = find(diff(SHOCK_frame_z_trace(:,1))<0); % 结束时间
    % 添加阴影
    px           = [t1, t2, t2, t1]; % 阴影的 x 坐标
    py           = [0,  0,  1,  1 ]; % 阴影的 y 坐标
    patch( px, py, [178 102 255]/255, 'FaceAlpha', 1, 'EdgeColor', 'none' ); % 添加阴影
    hold on
    for trialidx = 1:5

        plot(  mean(Activity_SHOCK_trial(:,trialidx,index_2_Regroup(17:27),mouseidx),3)*7,'LineWidth', 7, 'Color',Color_trial(trialidx,:)  )
        % pause
        hold on

    end
    axis off;
end


Activity_SHOCK_trial_MOs     = mean( reshape((  Activity_SHOCK_trial(:,:,[11 12],:) ), [75 50]) , 2);
Activity_SHOCK_trial_SSp_bfd = mean( reshape((  Activity_SHOCK_trial(:,:,[17 18],:) ), [75 50]) , 2);
Activity_SHOCK_trial_RSPagl  = mean( reshape((  Activity_SHOCK_trial(:,:,[13 14],:) ), [75 50]) , 2);
Activity_SHOCK_trial_VISpm   = mean( reshape((  Activity_SHOCK_trial(:,:,[15 16],:) ), [75 50]) , 2);
Activity_SHOCK_trial_RSpd    = mean( reshape((  Activity_SHOCK_trial(:,:,[5   6],:) ), [75 50]) , 2);

    figure( 615)
    t1           = find(diff(SHOCK_frame_z_trace(:,1))>0)+1; % 开始时间
    t2           = find(diff(SHOCK_frame_z_trace(:,1))<0); % 结束时间
    % 添加阴影
    px           = [t1, t2, t2, t1]; % 阴影的 x 坐标
    py           = [0.04,  0.04,  0.14,  0.14 ]; % 阴影的 y 坐标
    patch( px, py, [178 102 255]/255, 'FaceAlpha', 1, 'EdgeColor', 'none' ); % 添加阴影
    hold on
    for trialidx = 1:5

        plot(  Activity_SHOCK_trial_RSpd,'LineWidth', 7, 'Color','r' )
        % pause
        hold on
        plot(  mean(Activity_SHOCK,2),'LineWidth', 7, 'Color','k'  )
        % pause

    end
    axis off;



%------------------------------------------

figure(620)
for mouseidx=1:5
    for regionidx = 1:27
        plot(2*Activity_SHOCK_avg(:,regionidx,mouseidx)+0.2*regionidx,'LineWidth', 1, 'Color',color_inuse(regionidx,:))
        hold on
    end
end


figure(621)
Activity_SHOCK_avg_mean = mean(Activity_SHOCK_avg,3);
for regionidx = 1:27
    plot(2*Activity_SHOCK_avg_mean(:,regionidx)+0.2*regionidx,'LineWidth', 2, 'Color',color_inuse(regionidx,:))
    hold on
end

[M,I] = max(Activity_SHOCK_avg_mean)
[B,Idx] = sort(I)

figure(622)
for regionidx = 1:27
    plot(4*Activity_SHOCK_avg_mean(:,Idx(regionidx))+0.2*regionidx,'LineWidth', 2, 'Color',color_inuse(regionidx,:))
    hold on
    plot( SHOCK_frame_z_trace(:,1)*30,'k'  )
    hold on;
    text(80,0.21*regionidx+0.2, Cortex_Incommon_Sorted{Idx(regionidx)});
    hold on
end

figure(623)
for mouseidx=1:5
    for regionidx = 1:27
        plot(4*Activity_SHOCK_avg(:,Idx(regionidx),mouseidx)+0.2*regionidx,'LineWidth', 1, 'Color',color_inuse(regionidx,:))
        hold on
        plot( SHOCK_frame_z_trace(:,1)*30,'k'  )
        hold on;
        text(80,0.21*regionidx+0.2, Cortex_Incommon_Sorted{Idx(regionidx)});
        hold on
    end
end

for mouseidx=1

    Activity_SHOCK_avg_mean_mouse = diff(  Activity_SHOCK_avg(:,:,mouseidx)  );
    [M,I] = max(Activity_SHOCK_avg_mean_mouse)
    [B,Idx] = sort(I)

    figure(624)
    for regionidx = 1:27
        plot(20*(Activity_SHOCK_avg_mean_mouse(:,Idx(regionidx)))+0.2*regionidx,'LineWidth', 2, 'Color',color_inuse(regionidx,:))
        hold on
        plot( SHOCK_frame_z_trace(:,1)*30,'k'  )
        hold on;
        text(24,0.21*regionidx, [Cortex_Incommon_Sorted{Idx(regionidx)} ':' num2str(regionidx) '/' num2str(find(ismember(Mice_STRUCT_T(mouseidx).REGION_NAME, Cortex_Incommon_Sorted{Idx(regionidx)})))]);
        hold on
    end

    hold off

    figure(625)
    imshow(Mice_STRUCT_T(mouseidx).ATLAS90)

%     pause

end

%%=======================SHOCK=============================================
%%=======================SHOCK=============================================

%%========================CUE==============================================
%%========================CUE==============================================

close all
Activity_CUE_avg          = [];
Trace_CUE_avg             = [];
Zeropad_left              = 5;%10;
Zeropad_right             = 40;%50;

% CUE_MARKER                = zeros(245,1);
% CUE_MARKER(6:205,1)       = 0.15;


for mouseidx              = 1:5

    for regionidx         = 1:27


        Activity_CUE      = NaN( 200+Zeropad_left+Zeropad_right , 5 );
        Trace_CUE         = NaN(200+Zeropad_left+Zeropad_right,5);
        CUE_frame_z_trace = zeros(200+Zeropad_left+Zeropad_right,5);
        SPEED_unify_CUE   = NaN(200+Zeropad_left+Zeropad_right,5);

        for trialidx      = 1:5

            CUE_frame     = find(Mice_STRUCT_T(mouseidx).NEURON.CUE(logical(Frame_Inuse_Mouse{mouseidx}))==trialidx);

            CUE_frame_z   = [CUE_frame(1)-Zeropad_left:CUE_frame(1)-1, CUE_frame, CUE_frame(end)+1:CUE_frame(end)+Zeropad_right];

            temp          = Mice_STRUCT_T(mouseidx).NEURON.CUE(logical(Frame_Inuse_Mouse{mouseidx}));

            CUE_frame_z_trace(1:length(CUE_frame_z),trialidx) = .15*(temp(CUE_frame_z)>0);

            Activity_CUE(1:length(CUE_frame_z),trialidx) = single(Trace_Region_Avg(regionidx , CUE_frame_z , mouseidx))';

            Trace_CUE(1:length(CUE_frame_z),trialidx)    = single(Trace_unify_AVG(  mouseidx , CUE_frame_z  ))';

            SPEED_unify_CUE(1:length(CUE_frame_z),trialidx) = SPEED_unify(mouseidx,CUE_frame_z);

            %         pause

        end

        Activity_CUE_avg(:,regionidx,mouseidx) = nanmean(Activity_CUE,2);

        Activity_CUE_Region_Trial_Mouse(regionidx,:,:,mouseidx) = Activity_CUE;

        SPEED_unify_CUE_Trial_Mouse(:,:,mouseidx) = SPEED_unify_CUE;

        Trace_CUE_avg(:,:,mouseidx) = Trace_CUE;

    end


end




close all
for mouseidx=1:5

    figure(700+mouseidx)

    t1 = find(diff(CUE_frame_z_trace(:,mouseidx))>0)+1; % 开始时间
    t2 = find(diff(CUE_frame_z_trace(:,mouseidx))<0); % 结束时间
    t2=50
    % 添加阴影
    px = [t1, t2, t2, t1]; % 阴影的 x 坐标
    py = [0.05, 0.05, 0.25, 0.25]; % 阴影的 y 坐标
    patch(px, py, [204 229 255]/255, 'FaceAlpha', 1, 'EdgeColor', 'none'); % 添加阴影
    hold on
    for trialidx=1:5

        plot(Trace_CUE_avg(1:50,trialidx,mouseidx)*1+trialidx*0.025,'LineWidth', 3, 'Color',Color_trial(trialidx,:))
        % pause
        hold on

    end

    %     plot( SHOCK_frame_z_trace(:,1)*5,'k'  )

    axis off;
    % pause

end

Trace_CUE_avg_T2 = squeeze(  Trace_CUE_avg(:,2,:)  );
Trace_CUE_avg_T3 = squeeze(  Trace_CUE_avg(:,3,:)  );

close all
for mouseidx=1:5

    figure(710+mouseidx)

    plot( CUE_frame_z_trace(:,1)*5,'k' ,'LineStyle', '-.', 'LineWidth', 3 )
    hold on

    for trialidx=1:5

        plot(Trace_CUE_avg(:,trialidx,mouseidx)*7,'LineWidth', 3, 'Color',Color_trial(trialidx,:))

        hold on

    end

    axis off;
    % pause

end



close all



Activity_CUE_trial_MOTION_MO    = squeeze(  mean(Activity_CUE_Region_Trial_Mouse(index_2_Regroup(1:4  ),:,:,:) , 1)  );
Activity_CUE_trial_MOTION_SSP   = squeeze(  mean(Activity_CUE_Region_Trial_Mouse(index_2_Regroup(5:14 ),:,:,:) , 1)  );
Activity_CUE_trial_MOTION_RSP   = squeeze(  mean(Activity_CUE_Region_Trial_Mouse(index_2_Regroup(15:16),:,:,:) , 1)  );
Activity_CUE_trial_MOTION_VIS   = squeeze(  mean(Activity_CUE_Region_Trial_Mouse(index_2_Regroup(21:26),:,:,:) , 1)  );
Activity_CUE_trial_MOTION_VIS_w = squeeze(  mean(Activity_CUE_Region_Trial_Mouse(index_2_Regroup(19:27),:,:,:) , 1)  );

Activity_CUE_trial_MOTION_all   = squeeze(  mean(Activity_CUE_Region_Trial_Mouse(:,:,:,:) , 1)  );

Activity_CUE_trial_MOTION_MO_t4    = squeeze(  mean(Activity_CUE_Region_Trial_Mouse(index_2_Regroup(1:4  ),:,4,:) , 1)  );
Activity_CUE_trial_MOTION_SSP_t4   = squeeze(  mean(Activity_CUE_Region_Trial_Mouse(index_2_Regroup(5:14 ),:,4,:) , 1)  );
Activity_CUE_trial_MOTION_RSP_t4   = squeeze(  mean(Activity_CUE_Region_Trial_Mouse(index_2_Regroup(15:16),:,4,:) , 1)  );
Activity_CUE_trial_MOTION_VIS_t4   = squeeze(  mean(Activity_CUE_Region_Trial_Mouse(index_2_Regroup(21:26),:,4,:) , 1)  );
Activity_CUE_trial_MOTION_VIS_w_t4 = squeeze(  mean(Activity_CUE_Region_Trial_Mouse(index_2_Regroup(19:27),:,4,:) , 1)  );
Activity_CUE_trial_MOTION_all_t4   = squeeze(  mean(Activity_CUE_Region_Trial_Mouse(:,:,4,:) , 1)  );


Activity_CUE_trial_MOTION_MO_t5    = squeeze(  mean(Activity_CUE_Region_Trial_Mouse(index_2_Regroup(1:4  ),:,5,:) , 1)  );
Activity_CUE_trial_MOTION_SSP_t5   = squeeze(  mean(Activity_CUE_Region_Trial_Mouse(index_2_Regroup(5:14 ),:,5,:) , 1)  );
Activity_CUE_trial_MOTION_RSP_t5   = squeeze(  mean(Activity_CUE_Region_Trial_Mouse(index_2_Regroup(15:16),:,5,:) , 1)  );
Activity_CUE_trial_MOTION_VIS_t5   = squeeze(  mean(Activity_CUE_Region_Trial_Mouse(index_2_Regroup(21:26),:,5,:) , 1)  );
Activity_CUE_trial_MOTION_VIS_w_t5 = squeeze(  mean(Activity_CUE_Region_Trial_Mouse(index_2_Regroup(19:27),:,5,:) , 1)  );

Activity_CUE_trial_MOTION_all_t5   = squeeze(  mean(Activity_CUE_Region_Trial_Mouse(:,:,5,:) , 1)  );

for mouseidx=1:5

    fig  = figure(900+mouseidx)
    
    Activity_CUE_trial_MOTION_avg  = squeeze(  Activity_CUE_trial_MOTION_all(:,:,mouseidx) );
    Activity_CUE_Region_Trial_reg    = squeeze(  Activity_CUE_trial_MOTION_VIS_w(:,:,mouseidx) );

    plot( CUE_frame_z_trace(:,1)*5,'k' ,'LineStyle', '-.', 'LineWidth', 3 )
    hold on

    for trialidx=5
        

        plot(Activity_CUE_trial_MOTION_avg(:,trialidx)*7,'LineWidth', 2, 'Color',Color_trial(trialidx,:))
        hold on
%         plot(Trace_CUE_avg(:,trialidx,mouseidx)*7,'LineWidth', 3, 'Color','k')
        plot(  Activity_CUE_Region_Trial_reg(:,trialidx)*7,'LineWidth', 3, 'Color','k');
         
        hold on

     
    end

end


axis off;

hold off

regionidx
pause 

% end




close all

regiongrp(1) = 0
regiongrp(2) = 4
regiongrp(3) = 14
regiongrp(4) = 18
regiongrp(5) = 27

for regiongrpidx=1:4

    fig=figure(950+regiongrpidx)

    Activity_CUE_MOTION_Trial_MICE = squeeze(   Activity_CUE_Region_Trial_Mouse(index_2_Regroup(regiongrp(regiongrpidx)+1:regiongrp(regiongrpidx+1)),:,:,1:5)  );

%     Activity_CUE_MOTION_Trial_MICE = squeeze(   Activity_CUE_Region_Trial_Mouse(index_2_Regroup(21:26),:,:,1:5)  );

    Activity_CUE_MOTION_Trial_MICE_avg = squeeze( mean(Activity_CUE_MOTION_Trial_MICE,1) );
    Activity_CUE_MOTION_Trial_MICE_avg_avgmice = squeeze( mean(Activity_CUE_MOTION_Trial_MICE_avg,3) );

    plot( CUE_frame_z_trace(:,1)*5,'k' ,'LineStyle', '-.', 'LineWidth', 3 )
    hold on

    for trialidx=1:5
        
        plot(Activity_CUE_MOTION_Trial_MICE_avg_avgmice(:,trialidx)*7,'LineWidth', 2, 'Color',Color_trial(trialidx,:))
        hold on
     
    end

axis off;

end


close all

for regionidx=1:27

    fig=figure(800+regionidx)

    for mouseidx=1:5

        Activity_CUE_Region_Trial_M = squeeze(Activity_CUE_Region_Trial_Mouse(:,:,:,mouseidx));

        for trialidx=1:5

            subplot(5,5,trialidx+(mouseidx-1)*5)
            t1 = find(diff(CUE_frame_z_trace(:,mouseidx))>0)+1; % 开始时间
            t2 = find(diff(CUE_frame_z_trace(:,mouseidx))<0); % 结束时间
            % 添加阴影
            px = [t1, t2, t2, t1]; % 阴影的 x 坐标
            py = [0, 0, 1, 1]; % 阴影的 y 坐标
            patch(px, py, [204 229 255]/255, 'FaceAlpha', 1, 'EdgeColor', 'none'); % 添加阴影
            hold on
            plot(Trace_CUE_avg(:,trialidx,mouseidx)*7,'k.')
            hold on
            plot(Activity_CUE_Region_Trial_M(regionidx,:,trialidx)*7,'LineWidth', 2, 'Color',Color_trial(trialidx,:))
            hold on
            plot(SPEED_unify_CUE_Trial_Mouse(:,trialidx,mouseidx)*0.5,'blue-')
            hold on
            text(250,0.5, Cortex_Incommon_Sorted{regionidx});

        end

    end
    %     pause

end
axis off;




close all

regions = [ 18 ]%18%[3 4 10 11 17 18 23 24]

for mouseidx=1:5

    fig=figure(6000+mouseidx)

    %         fig=figure(6000+regionidx)

    %           subplot(1,5,mouseidx)

    Activity_CUE_Region_Trial_M = squeeze(Activity_CUE_Region_Trial_Mouse(:,:,:,mouseidx));

    for trialidx=4:5

        subplot(1,2,trialidx-3)
        t1 = find(diff(CUE_frame_z_trace(:,mouseidx))>0)+1; % 开始时间
        t2 = find(diff(CUE_frame_z_trace(:,mouseidx))<0); % 结束时间
        % 添加阴影
        px = [t1, t2, t2, t1]; % 阴影的 x 坐标
        py = [0, 0, 1.5, 1.5]; % 阴影的 y 坐标
        patch(px, py, [204 229 255]/255, 'FaceAlpha', 1, 'EdgeColor', 'none'); % 添加阴影
        hold on

        for regionidx=1:length(regions)

            %             plot( CUE_frame_z_trace(:,1)*5,'k' ,'LineStyle', '-.', 'LineWidth', 3 )
            %             hold on


            Trace_CUE_avg_diff = (Trace_CUE_avg(:,trialidx,mouseidx)*7);
            Activity_CUE_Region_Trial_M_diff = ( Activity_CUE_Region_Trial_M(regions(regionidx),:,trialidx)*7);


            plot(Trace_CUE_avg_diff,'LineStyle', '-.','LineWidth', 1, 'Color','k' )
            hold on
            plot(Activity_CUE_Region_Trial_M_diff,'LineWidth', 3, 'Color',color_inuse(regions(regionidx),:))
            hold on
            %         plot(SPEED_unify_CUE_Trial_Mouse(:,trialidx,mouseidx)*0.5+(trialidx-4)*0.5+0.5,'red-')
            hold on
            %             text(65,0.5, Cortex_Incommon_Sorted{regions(regionidx)});
            %             axis off
            %             grid off

            temp = Activity_CUE_Region_Trial_M(regions(regionidx),t1,:);
            %         temp_min = min(temp)
            %         temp_max = max(temp);
            %             xlim([0 60])
            %                     ylim([0.5-0.1 0.5+1])

        end

    end
    %     pause

end
axis off;




close all

regions = [18 ]%18%[3 4 10 11 17 18 23 24]

    for mouseidx=1:5

  fig=figure(6000+mouseidx)

    %         fig=figure(6000+regionidx)

              %           subplot(1,5,mouseidx)

        Activity_CUE_Region_Trial_M = squeeze(Activity_CUE_Region_Trial_Mouse(:,:,:,mouseidx));
      
%         subplot(1,2,trialidx-3)
                    t1 = find(diff(CUE_frame_z_trace(:,mouseidx))>0)+1; % 开始时间
            t2 = find(diff(CUE_frame_z_trace(:,mouseidx))<0); % 结束时间
            % 添加阴影
            px = [t1, t2, t2, t1]; % 阴影的 x 坐标
            py = [0, 0, 3, 3]; % 阴影的 y 坐标
            patch(px, py, [204 229 255]/255, 'FaceAlpha', 1, 'EdgeColor', 'none'); % 添加阴影
            hold on  
        for trialidx=4:5



     
for regionidx=1:length(regions)
%             plot( CUE_frame_z_trace(:,1)*5,'k' ,'LineStyle', '-.', 'LineWidth', 3 )
%             hold on


            Trace_CUE_avg_diff = (Trace_CUE_avg(:,trialidx,mouseidx)*7);
            Activity_CUE_Region_Trial_M_diff = ( Activity_CUE_Region_Trial_M(regions(regionidx),:,trialidx)*7);


            plot(Trace_CUE_avg_diff+(trialidx-4)*1,'LineStyle', '-.','LineWidth', 1, 'Color','k' )
            hold on
            plot(Activity_CUE_Region_Trial_M_diff+(trialidx-4)*1,'LineWidth', 3, 'Color',Color_trial(trialidx,:))
            hold on
            %         plot(SPEED_unify_CUE_Trial_Mouse(:,trialidx,mouseidx)*0.5+(trialidx-4)*0.5+0.5,'red-')
            hold on
%             text(65,0.5, Cortex_Incommon_Sorted{regions(regionidx)});
            %             axis off
            %             grid off

            temp = Activity_CUE_Region_Trial_M(regions(regionidx),t1,:);
            %         temp_min = min(temp)
            %         temp_max = max(temp);
%             xlim([0 60])
%                     ylim([0.5-0.1 0.5+1])

        end

    end
    %     pause

end
axis off;




close all
regions = [3 4 10 11 17 18 23 24]




     

      for mouseidx=1:5

    fig=figure(7000+mouseidx)     

        t1 = find(diff(CUE_frame_z_trace(:,mouseidx))>0)+1; % 开始时间
        t2 = find(diff(CUE_frame_z_trace(:,mouseidx))<0); % 结束时间
        % 添加阴影
        px = [t1, t2, t2, t1]; % 阴影的 x 坐标
        py = [0, 0, 1.5, 1.5]; % 阴影的 y 坐标
        patch(px, py, [204 229 255]/255, 'FaceAlpha', 1, 'EdgeColor', 'none'); % 添加阴影
        hold on

for regionidx=1:length(regions)

        Activity_CUE_Region_Trial_M = squeeze(Activity_CUE_Region_Trial_Mouse(:,:,:,mouseidx));

 for trialidx=5

        
       
        plot(Activity_CUE_Region_Trial_M(regions(regionidx),:,trialidx)*7+regionidx*0.1,'LineWidth', 3, 'Color',color_inuse(regions(regionidx),:))
        hold on
       
        text(65,0.5, Cortex_Incommon_Sorted{regions(regionidx)});
        axis off
        grid off

        temp = Activity_CUE_Region_Trial_M(regions(regionidx),t1,:);
%         temp_min = min(temp)
%         temp_max = max(temp);
%         xlim([0 60])
%         ylim([0.2 4])

    end

    end
%     pause

end
axis off;

%%========================CUE==============================================
%%========================CUE==============================================


%%========================Permutation Test SPEED===========================
%%========================Permutation Test SPEED===========================

%%===Analyze speed related neurons of M94

Trace_MotionNeuron_NegSpeed = {};
Trace_VisualNeuron_NegSpeed = {};
Trace_MotionNeuron_Speed    = {};
Trace_VisualNeuron_Speed    = {};

for mouseidx                      = 1:5


    Speed_frames                  = ~isnan(Mice_STRUCT_T(mouseidx).NEURON.Speed);
    Speed_Values                  = Mice_STRUCT_T(mouseidx).NEURON.Speed(Speed_frames);
    Trace4Speed                   = Mice_STRUCT_T(mouseidx).NEURON.Trace_unify(:,Speed_frames);
    Timediff                      = [0;diff(Mice_STRUCT_T(mouseidx).NEURON.TIME)'];

    Timediff4Speed                = Timediff(Speed_frames);
    
    Speed_frames_1                = Speed_Values>(  max(Speed_Values)/10  );

    figure(599)
    plot( Speed_frames_1)
    hold on;
    plot(Trace4Speed(1,:))
      
    ns                            = size(Trace4Speed,1);
    numShuffles                   = 500;
    numShuffles                   = numShuffles+1;
    
    NoFrames                      = size(Trace4Speed,2);

    shift                         = floor( rand(numShuffles,1) *NoFrames);
    shift(1,1)                    = 0;
    indexp                        = 1:NoFrames;
    
    T_spk_diff                    = Timediff4Speed;    
    
    T_spk_diff_Speed              = T_spk_diff( Speed_frames_1  );
    T_spk_diff_Other              = T_spk_diff( ~Speed_frames_1 );
    
    T_spk_diff_speed_sum          = sum(T_spk_diff_Speed);
    T_spk_diff_Other_sum          = sum(T_spk_diff_Other);

    FR_other                      = []
    FR_speed                      = []

    parfor idxnumShf              = 1:numShuffles
    
        indexs2                   = mod(indexp+shift(idxnumShf), NoFrames);
        indexs2(indexs2==0)       = NoFrames;
        Trace4Speed_shf           = Trace4Speed(:,indexs2);
    
        Trace4Speed_shf_speed     = Trace4Speed_shf(:,Speed_frames_1);
     
        Trace4Speed_shf_speed_sum = sum(Trace4Speed_shf_speed,2);     
    
        FR_speed(idxnumShf,:)     = Trace4Speed_shf_speed_sum/T_spk_diff_speed_sum;


        Trace4Speed_shf_other     = Trace4Speed_shf(:,~Speed_frames_1);
     
        Trace4Speed_shf_other_sum = sum(Trace4Speed_shf_other,2);     
    
        FR_other(idxnumShf,:)     = Trace4Speed_shf_other_sum/T_spk_diff_Other_sum;
          
        idxnumShf
    
    end
    
    FR_speed_other                = FR_speed - FR_other;
    FR_speed_other_Real           = FR_speed_other(1,:);
    infoP_speed                   = sum(repmat(FR_speed_other_Real,numShuffles-1,1) >= FR_speed_other(2:end,:),1)./(numShuffles-1);
       
    Num_speed_095                 = sum(infoP_speed>=0.95,1)
    Num_speed_005                 = sum(infoP_speed<=0.05,1)
    
    tt=(infoP_speed>=0.95)
    
    infoP_speed                   = infoP_speed';
    
    %%-------------------------------------------------------------------------
  [infoP_speed_Value,infoP_speed_Idx] = sort(infoP_speed,'descend');
    
  Neuron_idx4speed                 = infoP_speed>=0.95;    
  Neuron_idx4other                 = infoP_speed<=0.05;

  figure(600)

  imshow(uint8(255* ones(size(Mice_STRUCT_T(mouseidx).ATLAS2COMET,1),size(Mice_STRUCT_T(mouseidx).ATLAS2COMET,2)) ));
  hold on;

for regionidx=1:length(Mice_STRUCT_T(mouseidx).REGION_ID)

    BW = imbinarize(double(Mice_STRUCT_T(mouseidx).ATLAS2COMET ==regionidx));
    boundaries = bwboundaries(BW);


    if (~isempty(boundaries))
        for k = 1:1

            b = boundaries{k};
            plot(b(:,2), b(:,1), 'Color' , [0 0 0], 'LineWidth', 1); % Plot boundary with green line
            %     fill(b(:,2)+1, b(:,1)+1, Mice_STRUCT_T(1).COLOR(11,:),'FaceAlpha',0.5,'LineStyle',':');
            hold on

        end
        hold on
    end
    hold on
    regionidx
end

    hold on
    scatter( Mice_STRUCT_T(mouseidx).COORDINATES(Neuron_idx4speed,1) , Mice_STRUCT_T(mouseidx).COORDINATES(Neuron_idx4speed,2),8,[1 0 0] ,'filled' )
    hold on

    scatter( Mice_STRUCT_T(mouseidx).COORDINATES(Neuron_idx4other,1) , Mice_STRUCT_T(mouseidx).COORDINATES(Neuron_idx4other,2),8,[0 0 1] ,'filled' )
    hold on
    
figure(601)
plot(  mean( 50*Mice_STRUCT_T(mouseidx).NEURON.Trace_unify(Neuron_idx4speed,:),1)  )
hold on
plot(  Mice_STRUCT_T(mouseidx).NEURON.Speed )

figure(602)
plot(  mean( 50*Mice_STRUCT_T(mouseidx).NEURON.Trace_unify(Neuron_idx4other,:),1)  )
hold on
plot(  Mice_STRUCT_T(mouseidx).NEURON.Speed )


figure(603)

% region_motion_idx           = find(strcmp(Mice_STRUCT_T(mouseidx).REGION_NAME,'MOp1_R'))%[9]%[9 10 31 32];
% region_motion_neuron_other  = find(Neuron_idx4other&ismember(Mice_STRUCT_T(mouseidx).Neuron_REGION_ID, region_motion_idx));
% region_visual_idx           = find(strcmp(Mice_STRUCT_T(mouseidx).REGION_NAME,'VISp1_R'))%[29]%[29 30]
% region_visual_neuron_other  = find(Neuron_idx4other&ismember(Mice_STRUCT_T(mouseidx).Neuron_REGION_ID, region_visual_idx));

[~, index_MO]               = ismember(string(Mice_STRUCT_T(mouseidx).REGION_NAME) , {'MOp1_R','MOp1_L','MOs1_R','MOs1_L'});

region_motion_idx           = find(  index_MO )%[9]%[9 10 31 32];

region_motion_neuron_other  = find(Neuron_idx4other&ismember(Mice_STRUCT_T(mouseidx).Neuron_REGION_ID, region_motion_idx));

[~, index_VIS]              = ismember(string(Mice_STRUCT_T(mouseidx).REGION_NAME) , {'VISp1_R','VISp1_L','VISpm1_R','VISpm1_L','VISam1_R','VISam1_L'});

region_visual_idx           = find(  index_VIS )%[29]%[29 30]

region_visual_neuron_other  = find(Neuron_idx4other&ismember(Mice_STRUCT_T(mouseidx).Neuron_REGION_ID, region_visual_idx));


Trace_MotionNeuron_NegSpeed{mouseidx} = mean(Mice_STRUCT_T(mouseidx).NEURON.Trace_unify(region_motion_neuron_other,:) ,1 );
Trace_VisualNeuron_NegSpeed{mouseidx} = mean(Mice_STRUCT_T(mouseidx).NEURON.Trace_unify(region_visual_neuron_other,:) ,1 );

% Trace_MotionNeuron_NegSpeed = normalize(Trace_MotionNeuron_NegSpeed);
% Trace_VisualNeuron_NegSpeed = normalize(Trace_VisualNeuron_NegSpeed);

plot(   2*Trace_MotionNeuron_NegSpeed{mouseidx} ,'b'  )
hold on
plot(   2*Trace_VisualNeuron_NegSpeed{mouseidx}  ,'r' )
hold on
plot(Mice_STRUCT_T(mouseidx).NEURON.Speed/50 ,'k')


figure(604)

region_motion_neuron_speed = find(Neuron_idx4speed&ismember(Mice_STRUCT_T(mouseidx).Neuron_REGION_ID, region_motion_idx));

region_visual_neuron_speed = find(Neuron_idx4speed&ismember(Mice_STRUCT_T(mouseidx).Neuron_REGION_ID, region_visual_idx));

Trace_MotionNeuron_Speed{mouseidx} = mean(Mice_STRUCT_T(mouseidx).NEURON.Trace_unify(region_motion_neuron_speed,:) ,1 );
Trace_VisualNeuron_Speed{mouseidx} = mean(Mice_STRUCT_T(mouseidx).NEURON.Trace_unify(region_visual_neuron_speed,:) ,1 );

% Trace_MotionNeuron_Speed = normalize(Trace_MotionNeuron_Speed);
% Trace_VisualNeuron_Speed = normalize(Trace_VisualNeuron_Speed);

plot(   2*Trace_MotionNeuron_Speed{mouseidx} ,'b'  )
hold on
plot(   2*Trace_VisualNeuron_Speed{mouseidx} ,'r' )
hold on
plot(Mice_STRUCT_T(mouseidx).NEURON.Speed/50 ,'k')


%-------------How many neurons in each region are positively related to
%speed and how many negatively related


Percentage_Motion4speed(mouseidx) = sum(Neuron_idx4speed&ismember(Mice_STRUCT_T(mouseidx).Neuron_REGION_ID, region_motion_idx))/sum(ismember(Mice_STRUCT_T(mouseidx).Neuron_REGION_ID, region_motion_idx))
Percentage_Visual4speed(mouseidx) = sum(Neuron_idx4speed&ismember(Mice_STRUCT_T(mouseidx).Neuron_REGION_ID, region_visual_idx))/sum(ismember(Mice_STRUCT_T(mouseidx).Neuron_REGION_ID, region_visual_idx))
Percentage_Motion4other(mouseidx) = sum(Neuron_idx4other&ismember(Mice_STRUCT_T(mouseidx).Neuron_REGION_ID, region_motion_idx))/sum(ismember(Mice_STRUCT_T(mouseidx).Neuron_REGION_ID, region_motion_idx))
Percentage_Visual4other(mouseidx) = sum(Neuron_idx4other&ismember(Mice_STRUCT_T(mouseidx).Neuron_REGION_ID, region_visual_idx))/sum(ismember(Mice_STRUCT_T(mouseidx).Neuron_REGION_ID, region_visual_idx))

disp(['The mouse index: ' num2str(mouseidx)])

end


Percentage_Motion4speed = Percentage_Motion4speed'
Percentage_Visual4speed = Percentage_Visual4speed'
Percentage_Motion4other = Percentage_Motion4other'
Percentage_Visual4other = Percentage_Visual4other'

figure(605)

for mouseidx=1

    Trace_MotionNeuron_NegSpeed_mouse = Trace_MotionNeuron_NegSpeed{mouseidx};

    Trace_VisualNeuron_NegSpeed_mouse = Trace_VisualNeuron_NegSpeed{mouseidx};

    Trace_MotionNeuron_Speed_mouse    = Trace_MotionNeuron_Speed{mouseidx};
    Trace_VisualNeuron_Speed_mouse    = Trace_VisualNeuron_Speed{mouseidx};


    Trace_MotionNeuron_NegSpeed_mouse_FrameinUse(:,mouseidx) = Trace_MotionNeuron_NegSpeed_mouse( logical(Frame_Inuse_Mouse{1,mouseidx}) );

    Trace_VisualNeuron_NegSpeed_mouse_FrameinUse(:,mouseidx) = Trace_VisualNeuron_NegSpeed_mouse( logical(Frame_Inuse_Mouse{1,mouseidx}) );
    Trace_MotionNeuron_Speed_mouse_FrameinUse(:,mouseidx)    = Trace_MotionNeuron_Speed_mouse( logical(Frame_Inuse_Mouse{1,mouseidx}) );
    Trace_VisualNeuron_Speed_mouse_FrameinUse(:,mouseidx)    = Trace_VisualNeuron_Speed_mouse( logical(Frame_Inuse_Mouse{1,mouseidx}) );


    plot(   2*Trace_MotionNeuron_NegSpeed_mouse_FrameinUse(:,mouseidx) ,'b'  )
    hold on
    plot(   2*Trace_VisualNeuron_NegSpeed_mouse_FrameinUse(:,mouseidx) ,'r' )
    hold on

end

% plot(Mice_STRUCT_T(mouseidx).NEURON.Speed/50 ,'k')


figure(606)

for mouseidx=1
plot(   2*Trace_MotionNeuron_Speed_mouse_FrameinUse(:,mouseidx) ,'b'  )
hold on
plot(   2*Trace_VisualNeuron_Speed_mouse_FrameinUse(:,mouseidx) ,'r' )
hold on
end

% plot(Mice_STRUCT_T(mouseidx).NEURON.Speed/50 ,'k')
%-------------------------------------------------------------------------




%%========================Permutation Test SPEED===========================
%%========================Permutation Test SPEED===========================



%%========================PEARSON CORRELTION===============================
%%========================PEARSON CORRELTION===============================

close all

for mouseidx=1:5

   SHOCK_frame_T = find(Mice_STRUCT_T(mouseidx).NEURON.SHOCK(logical(Frame_Inuse_Mouse{mouseidx}))>0);

   idx=1;
   Activity_SHOCK_T=[]
    for regionidx = [3 4 5 6 11 12]
              
            Activity_SHOCK_T(:,idx) = single(Trace_Region_Avg(regionidx , SHOCK_frame_T , mouseidx))'; 
            idx = idx +1;
                            
    end

    corrMatrix = corr(Activity_SHOCK_T);

%     pause


figure(1100+mouseidx)
heatmap(corrMatrix);

end
%-------------START-Sliding window correlation-----------------------------




close all

    sliding_correlation_Trace   = []
    sliding_correlation_Cortex  = []


for mouseidx=1:5



for scale=12
window_size = 50*scale;  % 滑窗大小
step_size   = 1;     % 滑动步长
num_windows = floor((length(Cortex_unify_AVG) - window_size) / step_size) + 1;


% 计算每个滑窗的相关性
for i = 1:num_windows
    % 窗口的起始和结束索引
    start_idx = (i-1) * step_size + 1;
    end_idx = start_idx + window_size - 1;
    
    % 提取窗口内的信号
    window_signal1 = single( Trace_unify_AVG(mouseidx,start_idx:end_idx) )';
    window_signal3 = single( SPEED_unify(mouseidx,start_idx:end_idx) )';
    window_signal2 = single( Cortex_unify_AVG(mouseidx,start_idx:end_idx) )';
    % 计算窗口内两个信号的相关系数
    corr_matrix = corr(window_signal1, window_signal3, 'rows','complete');
    
    % 提取相关系数 (corrcoef 返回的是一个矩阵)
    sliding_correlation_Trace(i,mouseidx) = corr_matrix;

    corr_matrix = corr(window_signal2, window_signal3, 'rows','complete');
    
    % 提取相关系数 (corrcoef 返回的是一个矩阵)
    sliding_correlation_Cortex(i,mouseidx) = corr_matrix;

    i

end


FIG_CORRELATION = figure(500+mouseidx)
% 绘制相关系数随时间的变化
time_vector = (1:num_windows) * step_size;  % 计算滑窗位置对应的时间点
plot(time_vector, sliding_correlation_Cortex(:,mouseidx), '-', 'LineWidth', 3,'Color','k');
hold on;
plot(time_vector, sliding_correlation_Trace(:,mouseidx), '-', 'LineWidth', 3,'Color',[1 0.4 0.4]);
hold on;
line([1 7000] , [0 0 ] , 'Color','k','LineWidth', 1 , 'LineStyle','--')

xlabel('Time (Sliding Window Position)');
ylabel('Correlation Coefficient');
title('Correlation between Activity and Head Speed');
hold on;
CueShock_marker_color(miniscopecuebeginFrame, miniscopecueendFrame,miniscopeshockbeginFrame,miniscopeshockendFrame,0);

% plot(  (Mice_STRUCT_T(mouseidx).NEURON.SHOCK(logical(Frame_Inuse_Mouse{mouseidx}))>0)*2-1,'k'  )
% hold on;
% plot(  (Mice_STRUCT_T(mouseidx).NEURON.CUE(logical(Frame_Inuse_Mouse{mouseidx}))>0)*2-1,'k'    )
% hold on


FIG_Zscore = figure(520+mouseidx)

plot( normalize( SPEED_unify(mouseidx,:) )', '-','Color','k');
hold on;
plot( normalize( Trace_unify_AVG(mouseidx,:) )',  '-','Color','r');
hold on;
plot(  (Mice_STRUCT_T(mouseidx).NEURON.SHOCK(logical(Frame_Inuse_Mouse{mouseidx}))>0)*2-1,'k'  )
hold on;
plot(  (Mice_STRUCT_T(mouseidx).NEURON.CUE(logical(Frame_Inuse_Mouse{mouseidx}))>0)*2-1,'k'    )
hold on

FIG_Zscore = figure(540+mouseidx)

histogram( normalize( SPEED_unify(mouseidx,:) )');
hold on;
histogram( normalize( Trace_unify_AVG(mouseidx,:) )');
hold on;
histogram( normalize( Cortex_unify_AVG(mouseidx,:) )');



pause 

end


SHOCK_FRAME_1(:,mouseidx)  = double(Mice_STRUCT_T(mouseidx).NEURON.SHOCK(logical(Frame_Inuse_Mouse{mouseidx}))>0);
CUE_FRAME_1(  :,mouseidx)  = double(Mice_STRUCT_T(mouseidx).NEURON.CUE(  logical(Frame_Inuse_Mouse{mouseidx}))>0);

SHOCK_FRAME_n1(:,mouseidx) = -(Mice_STRUCT_T(mouseidx).NEURON.SHOCK(logical(Frame_Inuse_Mouse{mouseidx}))>0);
CUE_FRAME_n1(  :,mouseidx) = -(Mice_STRUCT_T(mouseidx).NEURON.CUE(  logical(Frame_Inuse_Mouse{mouseidx}))>0);

end


SHOCK_FRAME_1 = double( SHOCK_FRAME_1 );
CUE_FRAME_1   = double( CUE_FRAME_1   );


 %-------------START-Sliding window correlation-----------------------------



%-------------START-Cortex Regions Sliding window correlation-----------------------------


find(Mice_STRUCT_T(mouseidx).NEURON.SHOCK(logical(Frame_Inuse_Mouse{mouseidx}))>0)

Mice_STRUCT_T(mouseidx).NEURON.CUE(logical(Frame_Inuse_Mouse{mouseidx}))>0


%%==============FOR SLIDING CORRELATION OF MOTER REGIONS AND SPEED=========


close all

for mouseidx    = 4%1:5


    scale=12
    window_size = 50*scale;  % 滑窗大小
    step_size   = 1;     % 滑动步长
    num_windows = floor((length(Cortex_unify_AVG) - window_size) / step_size) + 1;

    % 初始化相关性存储数组
    sliding_correlation          = zeros(1, num_windows);

    Trace_unify_Avg_Matrix       = squeeze(  Trace_unify_AVG( mouseidx , : )  );

    Cortex_unify_Avg_Matrix      = squeeze( Cortex_unify_AVG( mouseidx , : )  );
    
    RawAct_unify_Avg_Matrix      = squeeze( RawAct_unify_AVG( mouseidx , : )  );

    sliding_correlation_Trace    = []
    sliding_correlation_TraceReg = []
    sliding_correlation_CortexReg = []




%----------------Motion regions

%     for regionidx                = [3,4,5,6,7,8,11,12, 13,14,17,18,19,20,21,22,23,24,25,26,27] ;%1:size(Cortex_Incommon,1)  %[3,4,7,21,22] for motion related and %[9:10,15:16] for visualiation related
%  for regionidx                = [3,4,11,12] ;
%----------------Motion regions


%----------------Visual regions

    for regionidx                = [1 2 9 10 15 16];%[1,2,9,10,15,16,25:27];%1:size(Cortex_Incommon,1)  %[3,4,7,21,22] for motion related and %[9:10,15:16] for visualiation related

%----------------Visual regions

        Cortex_Region_Avg_Matrix = squeeze(  Cortex_Region_Avg(regionidx , : , mouseidx)  );
        Trace_Region_Avg_Matrix  = squeeze(  Trace_Region_Avg( regionidx , : , mouseidx)  );
       
        RawAct_Region_Avg_Matrix = squeeze( RawAct_Region_Avg( regionidx , : , mouseidx)  );

        % 计算每个滑窗的相关性
        for i              = 1:num_windows
            % 窗口的起始和结束索引
            start_idx      = (i-1) * step_size + 1;
            end_idx        = start_idx + window_size - 1;

            % 提取窗口内的信号
            window_signal1 = single( Trace_Region_Avg_Matrix(start_idx:end_idx)   )';
            window_signal3 = single( SPEED_unify(mouseidx,start_idx:end_idx)      )';
            window_signal2 = single( Trace_unify_Avg_Matrix( start_idx:end_idx)   )';

            window_signal4 = single( Cortex_unify_Avg_Matrix( start_idx:end_idx)  )';
            window_signal5 = single( Cortex_Region_Avg_Matrix( start_idx:end_idx) )';

            window_signal6 = single( RawAct_Region_Avg_Matrix( start_idx:end_idx) )';
            window_signal7 = single( RawAct_unify_Avg_Matrix( start_idx:end_idx)  )';
            % 计算窗口内两个信号的相关系数
            corr_matrix    = corr(window_signal1, window_signal3, 'rows','complete');
            sliding_correlation_TraceReg(regionidx,i) = corr_matrix;

            corr_matrix    = corr(window_signal2, window_signal3, 'rows','complete');       
            sliding_correlation_TraceAll(regionidx,i) = corr_matrix;

            corr_matrix    = corr(window_signal4, window_signal3, 'rows','complete');
            sliding_correlation_CortexAll(regionidx,i) = corr_matrix;

            corr_matrix    = corr(window_signal5, window_signal3, 'rows','complete');       
            sliding_correlation_CortexReg(regionidx,i) = corr_matrix;

            corr_matrix    = corr(window_signal6, window_signal3, 'rows','complete');       
            sliding_correlation_CortexROIReg(regionidx,i) = corr_matrix;

            corr_matrix    = corr(window_signal7, window_signal3, 'rows','complete');       
            sliding_correlation_CortexROIAll(regionidx,i) = corr_matrix;

            i

        end


%         figure(600+10*mouseidx+1)
%         set(gcf, 'Position', [100, 100, 600, 380]);      
%         
%         CueShock_marker_color(miniscopecuebeginFrame, miniscopecueendFrame,miniscopeshockbeginFrame,miniscopeshockendFrame,0);

        % 绘制相关系数随时间的变化

           figure(600+10*mouseidx+1)
    set(gcf, 'Position', [100, 100, 600, 380]);      
        
    CueShock_marker_color(miniscopecuebeginFrame, miniscopecueendFrame,miniscopeshockbeginFrame,miniscopeshockendFrame,0);

        time_vector     = (1:num_windows) * step_size;  % 计算滑窗位置对应的时间点
        plot(time_vector, sliding_correlation_TraceAll(regionidx,:) , '-', 'LineWidth', 2.5, 'Color','k');
        hold on;
        plot(time_vector, sliding_correlation_TraceReg(regionidx,:)  , '-', 'LineWidth' ,1.5, 'Color',color_inuse(round(regionidx*1),:));%,'DisplayName',Cortex_Incommon_Sorted{regionidx}
        hold on;
        line([1 7000] , [0 0] , 'Color','k','LineWidth', 1 , 'LineStyle','--')
        box off
   
        xlim([0 7000])
        ylim([-1 1])
        xticks([0:1000:7000])
        yticks([-1:0.5:1])
        % hold on
        % text(6200,1.5*regionidx-1.5, Cortex_Incommon_Sorted{regionidx});
        set(gca, 'TickDir', 'out') 
        set(gca, 'LineWidth', 1.5) 
        xlabel('Time (Sliding Window Position)');
        ylabel('Correlation');
        title([Mice_STRUCT_T(mouseidx).NAME(4:end) ]);%': ' Cortex_Incommon_Sorted{regionidx}
        set(gca, 'FontSize',15)
        set(gca, 'XTick', []);
        set(gca, 'box', 'off', 'XColor', 'none');
        hold on
        
  


         figure(700+10*mouseidx+1)
    set(gcf, 'Position', [100, 100, 600, 380]);      
        
    CueShock_marker_color(miniscopecuebeginFrame, miniscopecueendFrame,miniscopeshockbeginFrame,miniscopeshockendFrame,0);

        time_vector     = (1:num_windows) * step_size;  % 计算滑窗位置对应的时间点
        plot(time_vector, sliding_correlation_CortexAll(regionidx,:) , '-', 'LineWidth', 2.5, 'Color','k');
        hold on;
        plot(time_vector, sliding_correlation_CortexReg(regionidx,:)  , '-', 'LineWidth' ,1.5, 'Color',color_inuse(round(regionidx*1),:));%,'DisplayName',Cortex_Incommon_Sorted{regionidx}
        hold on;
        line([1 7000] , [0 0] , 'Color','k','LineWidth', 1 , 'LineStyle','--')
        box off
   
        xlim([0 7000])
        ylim([-1 1])
        xticks([0:1000:7000])
        yticks([-1:0.5:1])
        % hold on
        % text(6200,1.5*regionidx-1.5, Cortex_Incommon_Sorted{regionidx});
        set(gca, 'TickDir', 'out') 
        set(gca, 'LineWidth', 1.5) 
        xlabel('Time (Sliding Window Position)');
        ylabel('Correlation');
        title([Mice_STRUCT_T(mouseidx).NAME(4:end) ]);%': ' Cortex_Incommon_Sorted{regionidx}
        set(gca, 'FontSize',15)
        set(gca, 'XTick', []);
        set(gca, 'box', 'off', 'XColor', 'none');
        hold on
        

      figure(800+10*mouseidx+1)
    set(gcf, 'Position', [100, 100, 600, 380]);      
        
    CueShock_marker_color(miniscopecuebeginFrame, miniscopecueendFrame,miniscopeshockbeginFrame,miniscopeshockendFrame,0);

        time_vector     = (1:num_windows) * step_size;  % 计算滑窗位置对应的时间点
        plot(time_vector, sliding_correlation_CortexROIAll(regionidx,:) , '-', 'LineWidth', 2.5, 'Color','k');
        hold on;
        plot(time_vector, sliding_correlation_CortexROIReg(regionidx,:)  , '-', 'LineWidth' ,1.5, 'Color',color_inuse(round(regionidx*1),:));%,'DisplayName',Cortex_Incommon_Sorted{regionidx}
        hold on;
        line([1 7000] , [0 0] , 'Color','k','LineWidth', 1 , 'LineStyle','--')
        box off
   
        xlim([0 7000])
        ylim([-1 1])
        xticks([0:1000:7000])
        yticks([-1:0.5:1])
        % hold on
        % text(6200,1.5*regionidx-1.5, Cortex_Incommon_Sorted{regionidx});
        set(gca, 'TickDir', 'out') 
        set(gca, 'LineWidth', 1.5) 
        xlabel('Time (Sliding Window Position)');
        ylabel('Correlation');
        title([Mice_STRUCT_T(mouseidx).NAME(4:end) ]);%': ' Cortex_Incommon_Sorted{regionidx}
        set(gca, 'FontSize',15)
        set(gca, 'XTick', []);
        set(gca, 'box', 'off', 'XColor', 'none');
        hold on
        
        %          plot(  (Mice_STRUCT_T(mouseidx).NEURON.SHOCK(logical(Frame_Inuse_Mouse{mouseidx}))>0)*2-1,'k'  )
        %          hold on;
        %          plot(  (Mice_STRUCT_T(mouseidx).NEURON.CUE(logical(Frame_Inuse_Mouse{mouseidx}))>0)*2-1,'k'    )
        %          hold on

%         pause


    end




%     hold off


end


sliding_correlation_Trace = sliding_correlation_Trace';
%%==============FOR SLIDING CORRELATION OF MOTER REGIONS AND SPEED=========
%-------------START-MOP Neurons Sliding window correlation-----------------



%--------------------------------------------------------------------------
%-------------START-Motion regions SPON VS TRAINING CORRELATION------------
%-------------Compare average of neurons and cortex to speed

close all
    Correlation_CortexReg_T = []
    Correlation_TraceReg_T = []
Correlation_RawActReg_T=[]

for mouseidx    = 1:5

%     MotionRegionsID = [3,4,5,6,7,8,11,12, 13,14,17,18,19,20,21,22,23,24,25,26,27] %[3,4,7,8,11,12,21,22];
    MOID            = [3:4,11:12]
    MOp             = [3:4]
    MOs             = [11:12]
    SSPID           = [7:8,17:24]
    VISID           = [1:2,9:10,15:16]
    VISID2          = [25:27]
    RSPID           = [5:6,13:14]
    WholeRegionsID  = 1:27;
    RegionsChosen   = MOp;

    for regionidx                = 1:length(RegionsChosen) %1:size(Cortex_Incommon,1)  %[3,4,7,21,22] for motion related and %[9:10,15:16] for visualiation related


        Trace_Region_Avg_Matrix  = squeeze(  Trace_Region_Avg(  RegionsChosen(regionidx) , : , mouseidx)  );
        Cortex_Region_Avg_Matrix = squeeze(  Cortex_Region_Avg( RegionsChosen(regionidx) , : , mouseidx)  );
        RawAct_Region_Avg_Matrix = squeeze(  RawAct_Region_Avg( RegionsChosen(regionidx) , : , mouseidx)  );

        Trace_SPON     = single(  Trace_Region_Avg_Matrix( 1:miniscopecuebeginFrame(1)-1)  )';
        Cortex_SPON    = single(  Cortex_Region_Avg_Matrix(1:miniscopecuebeginFrame(1)-1)  )';
        RawAct_SPON    = single(  RawAct_Region_Avg_Matrix(1:miniscopecuebeginFrame(1)-1)  )';
        SPEED_SPON     = single(  SPEED_unify(mouseidx,    1:miniscopecuebeginFrame(1)-1)  )';

        Trace_TRAIN1   = single( Trace_Region_Avg_Matrix(  miniscopeshockendFrame(1):miniscopecuebeginFrame(2)-1)     )';
        Cortex_TRAIN1  = single( Cortex_Region_Avg_Matrix( miniscopeshockendFrame(1):miniscopecuebeginFrame(2)-1)     )';
        RawAct_TRAIN1  = single( RawAct_Region_Avg_Matrix( miniscopeshockendFrame(1):miniscopecuebeginFrame(2)-1)     )';
        SPEED_TRAIN1   = single( SPEED_unify(mouseidx,     miniscopeshockendFrame(1):miniscopecuebeginFrame(2)-1)     )';
        
        Trace_TRAIN2   = single( Trace_Region_Avg_Matrix(  miniscopeshockendFrame(2):miniscopecuebeginFrame(3)-1)     )';
        Cortex_TRAIN2  = single( Cortex_Region_Avg_Matrix( miniscopeshockendFrame(2):miniscopecuebeginFrame(3)-1)     )';
        RawAct_TRAIN2  = single( RawAct_Region_Avg_Matrix( miniscopeshockendFrame(2):miniscopecuebeginFrame(3)-1)     )';
        SPEED_TRAIN2   = single( SPEED_unify(mouseidx,     miniscopeshockendFrame(2):miniscopecuebeginFrame(3)-1)     )';
        
        Trace_TRAIN3   = single( Trace_Region_Avg_Matrix(  miniscopeshockendFrame(3):miniscopecuebeginFrame(4)-1)     )';
        Cortex_TRAIN3  = single( Cortex_Region_Avg_Matrix( miniscopeshockendFrame(3):miniscopecuebeginFrame(4)-1)     )';
        RawAct_TRAIN3  = single( RawAct_Region_Avg_Matrix( miniscopeshockendFrame(3):miniscopecuebeginFrame(4)-1)     )';
        SPEED_TRAIN3   = single( SPEED_unify(mouseidx,     miniscopeshockendFrame(3):miniscopecuebeginFrame(4)-1)     )';
        
        Trace_TRAIN4   = single( Trace_Region_Avg_Matrix(  miniscopeshockendFrame(4):miniscopecuebeginFrame(5)-1)     )';
        Cortex_TRAIN4  = single( Cortex_Region_Avg_Matrix( miniscopeshockendFrame(4):miniscopecuebeginFrame(5)-1)     )';
        RawAct_TRAIN4  = single( RawAct_Region_Avg_Matrix( miniscopeshockendFrame(4):miniscopecuebeginFrame(5)-1)     )';
        SPEED_TRAIN4   = single( SPEED_unify(mouseidx,     miniscopeshockendFrame(4):miniscopecuebeginFrame(5)-1)     )';

        Trace_TRAIN5   = single( Trace_Region_Avg_Matrix(  miniscopeshockendFrame(5):miniscopeshockendFrame(5)+538)     )';
        Cortex_TRAIN5  = single( Cortex_Region_Avg_Matrix( miniscopeshockendFrame(5):miniscopeshockendFrame(5)+538)     )';
        RawAct_TRAIN5  = single( RawAct_Region_Avg_Matrix( miniscopeshockendFrame(5):miniscopeshockendFrame(5)+538)     )';
        SPEED_TRAIN5   = single( SPEED_unify(mouseidx,     miniscopeshockendFrame(5):miniscopeshockendFrame(5)+538)     )';
        
        
        Correlation_TraceReg_T(1,regionidx,mouseidx)    = corr(SPEED_SPON, Trace_SPON,   'rows','complete');
        Correlation_TraceReg_T(2,regionidx,mouseidx)    = corr(SPEED_TRAIN1, Trace_TRAIN1, 'rows','complete');
        Correlation_TraceReg_T(3,regionidx,mouseidx)    = corr(SPEED_TRAIN2, Trace_TRAIN2, 'rows','complete');
        Correlation_TraceReg_T(4,regionidx,mouseidx)    = corr(SPEED_TRAIN3, Trace_TRAIN3, 'rows','complete');
        Correlation_TraceReg_T(5,regionidx,mouseidx)    = corr(SPEED_TRAIN4, Trace_TRAIN4, 'rows','complete');
        Correlation_TraceReg_T(6,regionidx,mouseidx)    = corr(SPEED_TRAIN5, Trace_TRAIN5, 'rows','complete');
       

        Correlation_CortexReg_T(1,regionidx,mouseidx)    = corr(SPEED_SPON, Cortex_SPON,   'rows','complete');
        Correlation_CortexReg_T(2,regionidx,mouseidx)    = corr(SPEED_TRAIN1, Cortex_TRAIN1, 'rows','complete');
        Correlation_CortexReg_T(3,regionidx,mouseidx)    = corr(SPEED_TRAIN2, Cortex_TRAIN2, 'rows','complete');
        Correlation_CortexReg_T(4,regionidx,mouseidx)    = corr(SPEED_TRAIN3, Cortex_TRAIN3, 'rows','complete');
        Correlation_CortexReg_T(5,regionidx,mouseidx)    = corr(SPEED_TRAIN4, Cortex_TRAIN4, 'rows','complete');
        Correlation_CortexReg_T(6,regionidx,mouseidx)    = corr(SPEED_TRAIN5, Cortex_TRAIN5, 'rows','complete');


        Correlation_RawActReg_T(1,regionidx,mouseidx)    = corr(SPEED_SPON,   RawAct_SPON,   'rows','complete');
        Correlation_RawActReg_T(2,regionidx,mouseidx)    = corr(SPEED_TRAIN1, RawAct_TRAIN1, 'rows','complete');
        Correlation_RawActReg_T(3,regionidx,mouseidx)    = corr(SPEED_TRAIN2, RawAct_TRAIN2, 'rows','complete');
        Correlation_RawActReg_T(4,regionidx,mouseidx)    = corr(SPEED_TRAIN3, RawAct_TRAIN3, 'rows','complete');
        Correlation_RawActReg_T(5,regionidx,mouseidx)    = corr(SPEED_TRAIN4, RawAct_TRAIN4, 'rows','complete');
        Correlation_RawActReg_T(6,regionidx,mouseidx)    = corr(SPEED_TRAIN5, RawAct_TRAIN5, 'rows','complete');
        


    end


end

Correlation_TraceReg_T_2D  = reshape(Correlation_TraceReg_T , [size(Correlation_TraceReg_T,1) size(Correlation_TraceReg_T,2)*size(Correlation_TraceReg_T,3)]);
Correlation_CortexReg_T_2D = reshape(Correlation_RawActReg_T , [size(Correlation_RawActReg_T,1) size(Correlation_RawActReg_T,2)*size(Correlation_RawActReg_T,3)]);


%-------------START-Motion regions SPON VS TRAINING CORRELATION------------
%--------------------------------------------------------------------------


close all

for mouseidx=1:5


    scale       = 12
    window_size = 50*scale;  % 滑窗大小
    step_size   = 10;     % 滑动步长
    num_windows = floor((length(Cortex_unify_AVG) - window_size) / step_size) + 1;

    % 初始化相关性存储数组
    sliding_correlation = zeros(1, num_windows);

    Trace_unify_Avg_Matrix   = squeeze(  Trace_unify_AVG( mouseidx , : )  );

    sliding_correlation_Trace=[]
    sliding_correlation_TraceReg=[]



    regionidx            = 3 %[9:10,15:16]%
    Trace_Neuron_ID(:,:) = Mice_STRUCT_T(mouseidx).NEURON.Trace_unify(Mice_STRUCT_T(mouseidx).Neuron_REGION_ID==CortexID_Mouse_Sorted(mouseidx,regionidx),logical(Frame_Inuse_Mouse{mouseidx}));

    Cortex_Region_Avg_Matrix = squeeze(  Cortex_Region_Avg(regionidx , : , mouseidx)  );
    Trace_Region_Avg_Matrix  = squeeze(  Trace_Region_Avg( regionidx , : , mouseidx)  );

    for neuronidx=1:size(Trace_Neuron_ID,1)

        Trace_Neuron_ID_Matrix = Trace_Neuron_ID(neuronidx,:) ;


        % 计算每个滑窗的相关性
        for i = 1:num_windows
            % 窗口的起始和结束索引
            start_idx = (i-1) * step_size + 1;
            end_idx = start_idx + window_size - 1;

            % 提取窗口内的信号
            window_signal2 = single( Trace_Region_Avg_Matrix(start_idx:end_idx) )';
            window_signal3 = single( SPEED_unify(mouseidx,start_idx:end_idx) )';
            window_signal1 = single( Trace_Neuron_ID_Matrix( start_idx:end_idx) )';
            % 计算窗口内两个信号的相关系数
            corr_matrix = corr(window_signal1, window_signal3, 'rows','complete');

            % 提取相关系数 (corrcoef 返回的是一个矩阵)
            sliding_correlation_TraceReg(i) = corr_matrix;

            corr_matrix = corr(window_signal2, window_signal3, 'rows','complete');

            % 提取相关系数 (corrcoef 返回的是一个矩阵)
            sliding_correlation_Trace(i) = corr_matrix;

            i

        end


        FIG_CORRELATION = figure(700+mouseidx)
        % 绘制相关系数随时间的变化
        time_vector = (1:num_windows) * step_size;  % 计算滑窗位置对应的时间点
        plot(time_vector, sliding_correlation_Trace , '-','Color','k');
        hold on;
        plot(time_vector, sliding_correlation_TraceReg  , '-');%,'DisplayName',Cortex_Incommon_Sorted{regionidx}


        xlabel('Time (Sliding Window Position)');
        ylabel('Correlation Coefficient');
        title(['Cortex Region: ' Cortex_Incommon_Sorted{regionidx} 'Correlation']);

        CueShock_marker_color(miniscopecuebeginFrame, miniscopecueendFrame,miniscopeshockbeginFrame,miniscopeshockendFrame,0);
        hold off
        
        neuronidx

        pause

    end

   

    hold off

    pause

end

%%========================PEARSON CORRELTION===============================
%%========================PEARSON CORRELTION===============================


%-------------START-MOP Neurons Sliding window correlation-----------------


for mouseidx = 1

    
    for regionidx                = 1:size(CortexID_Mouse_Sorted,2)

        Trace_Neuron_ID              = []

        Trace_Neuron_ID(:,:)     = Mice_STRUCT_T(mouseidx).NEURON.Trace_unify(Mice_STRUCT_T(mouseidx).Neuron_REGION_ID==CortexID_Mouse_Sorted(mouseidx,regionidx),logical(Frame_Inuse_Mouse{mouseidx}));

        Neuron_Cortex_Chosen_Num = Mice_STRUCT_T(mouseidx).REGION_ID( CortexID_Mouse_Sorted(mouseidx,regionidx) );

        Trace_Cortex_Chosen      = Mice_STRUCT_T(mouseidx).NEURON.Trace_unify_Region_Avg(CortexID_Mouse_Sorted(mouseidx,regionidx),logical(Frame_Inuse_Mouse{mouseidx}));

        Trace_Neuron_ID_Shocktrial_Sum = []

        for neuronidx = 1:Neuron_Cortex_Chosen_Num

            for trialidx = 1:5

                Shock_trial = (Mice_STRUCT_T(mouseidx).NEURON.SHOCK(logical(Frame_Inuse_Mouse{mouseidx}))==trialidx);

                Trace_Neuron_ID_Shocktrial_Sum(neuronidx,trialidx)  =  sum(  Trace_Neuron_ID(neuronidx,Shock_trial)  );

            end

        end


        for trialidx=1:5


            [Trace_Neuron_ID_Shock5th_Sum_sort , Index] = sort(Trace_Neuron_ID_Shocktrial_Sum(:,trialidx), 'descend');

            if (trialidx==1)

                Trace_Neuron_ID_Shock5th_Sum_sort_t1 = Trace_Neuron_ID_Shock5th_Sum_sort;
                
            end


            Trace_Neuron_ID_Shock5th_Chosen_Num(regionidx,trialidx) = sum(  Trace_Neuron_ID_Shock5th_Sum_sort  >  mean(Trace_Neuron_ID_Shock5th_Sum_sort_t1)  );

            Trace_Neuron_ID_Shock5th_Chosen_Intensity{regionidx} = Trace_Neuron_ID_Shock5th_Sum_sort;

            Trace_Neuron_ID_Shock5th_Sum_sort_per(regionidx,trialidx,mouseidx) = sum( Trace_Neuron_ID_Shock5th_Sum_sort  >  mean(Trace_Neuron_ID_Shock5th_Sum_sort_t1))/length(Trace_Neuron_ID_Shock5th_Sum_sort);


        end


        if strcmp(string(Cortex_Incommon_Sorted(regionidx)) , 'MOs1_L')

            figure(600+10*mouseidx)


            for neuronidx=1:Neuron_Cortex_Chosen_Num


                plot(  1*Trace_Neuron_ID(Index(neuronidx),:) + (neuronidx)-1 , 'LineWidth',1  )
                hold on;

                plot(  neuronidx*(Mice_STRUCT_T(mouseidx).NEURON.SHOCK(logical(Frame_Inuse_Mouse{mouseidx}))>0),'k'  )
                hold on;

                plot(  neuronidx*(Mice_STRUCT_T(mouseidx).NEURON.CUE(logical(Frame_Inuse_Mouse{mouseidx}))>0),'k'    )
                hold on;

                plot(  1* Trace_unify_SUM ./Trace_unify_NUM + 1*neuronidx-1,'g'  )
                hold on;


                plot( 1*Trace_Cortex_Chosen/Neuron_Cortex_Chosen_Num + 1*neuronidx-1,'b' )
                hold on

                plot( SPEED_unify + 1*neuronidx-1, 'LineWidth',1,'Color','k' )
                hold on
                text(6900,1*regionidx-1,strjoin([ string(Cortex_Incommon_Sorted(regionidx)) '=' num2str(Trace_Neuron_ID_Shock5th_Chosen_Num(regionidx,trialidx)) '/' num2str(Neuron_Cortex_Chosen_Num) 'neurons' ]));


                Trace_Neuron_ID_Shocktrial_Mean(regionidx,:) = mean( Trace_Neuron_ID( Index( Trace_Neuron_ID_Shock5th_Sum_sort  >  mean(Trace_Neuron_ID_Shock5th_Sum_sort) ),:) );

                chosen(regionidx,:) = mean( Trace_Neuron_ID_Shocktrial_Sum( Index( Trace_Neuron_ID_Shock5th_Sum_sort  >  mean(Trace_Neuron_ID_Shock5th_Sum_sort_t1) ),:) );


%                 pause

            end


        end


    end


end




figure(  700 + 10*mouseidx  )

chosen_avg = mean(chosen,1,'omitnan');
chosen_std = std(chosen,1,'omitnan');

x          = 1:size(chosen,2);

fill(  [x, flip(x)], [chosen_avg+chosen_std, flip(chosen_avg-chosen_std)], [0.8 0.8 0.8]  )
hold on;
plot(chosen_avg,'r')

title(Files(mouseidx).name)

figure(800+10*mouseidx)

plot(    Trace_Neuron_ID_Shock5th_Sum_sort_per(:,:,mouseidx)    )
xticks(1:27)
xticklabels(  C_sort(:)  )
title(Files(mouseidx).name)


figure(1000+10*mouseidx)

Trace_Neuron_ID_Shock5th_Chosen_Num_avg = mean(Trace_Neuron_ID_Shock5th_Chosen_Num,1,'omitnan');
Trace_Neuron_ID_Shock5th_Chosen_Num_std = std(Trace_Neuron_ID_Shock5th_Chosen_Num,1,'omitnan');
x                                       = 1:size(Trace_Neuron_ID_Shock5th_Chosen_Num,2);

fill(  [x, flip(x)], [Trace_Neuron_ID_Shock5th_Chosen_Num_avg+Trace_Neuron_ID_Shock5th_Chosen_Num_std, flip(Trace_Neuron_ID_Shock5th_Chosen_Num_avg-Trace_Neuron_ID_Shock5th_Chosen_Num_std)], [0.8 0.8 0.8]  )
hold on;

plot(Trace_Neuron_ID_Shock5th_Chosen_Num_avg)




figure(1200+10*mouseidx)

Trace_Neuron_ID_Shock5th_Chosen_Num_avg = mean(chosen ./ Trace_Neuron_ID_Shock5th_Chosen_Num,1,'omitnan');
Trace_Neuron_ID_Shock5th_Chosen_Num_std = std(chosen ./ Trace_Neuron_ID_Shock5th_Chosen_Num,1,'omitnan');
x = 1:size(Trace_Neuron_ID_Shock5th_Chosen_Num,2);

fill(  [x, flip(x)], [Trace_Neuron_ID_Shock5th_Chosen_Num_avg+Trace_Neuron_ID_Shock5th_Chosen_Num_std, flip(Trace_Neuron_ID_Shock5th_Chosen_Num_avg-Trace_Neuron_ID_Shock5th_Chosen_Num_std)], [0.8 0.8 0.8]  )
hold on;

plot(Trace_Neuron_ID_Shock5th_Chosen_Num_avg)




mouseidx


%          pause



%      end






for regionidx                = 1:size(bi_mouse_sorted,2)



    Trace_Neuron_ID_Shock5th_Sum_sort_per_region =  squeeze(  Trace_Neuron_ID_Shock5th_Sum_sort_per(regionidx,:,:)   );

    figure(1500+regionidx)

    plot(  Trace_Neuron_ID_Shock5th_Sum_sort_per_region  )
    title(  C_sort(regionidx)    )

    hold on


end


temp = squeeze(Trace_Neuron_ID_Shock5th_Sum_sort_per(:,:,4))




for mouseidx=1:5

    figure(1400+mouseidx)

    for regionidx =1:size(bi_mouse_sorted,2)

        Trace_Cortex_Chosen = Mice_STRUCT_T(mouseidx).NEURON.Trace_unify_Region_Avg(bi_mouse_sorted(mouseidx,regionidx),logical(Frame_Inuse_Mouse{mouseidx}));

        Neuron_Cortex_Chosen_Num = Mice_STRUCT_T(mouseidx).REGION_ID( bi_mouse_sorted(mouseidx,regionidx) );

        plot( 60*diff(Trace_Cortex_Chosen) + 1.5*regionidx-1.5 , 'LineWidth',1.5,'Color', color_inuse(regionidx,:))
        hold on;
        regionidx
        %             pause
        Trace_Cortex_Chosen = Trace_Cortex_Chosen*Neuron_Cortex_Chosen_Num;

        Neuron_Cortex_Chosen_Num_Sum(regionidx) = Neuron_Cortex_Chosen_Num_Sum(regionidx) + Neuron_Cortex_Chosen_Num;
        Trace_Cortex_Chosen_Sum(regionidx,:) = Trace_Cortex_Chosen_Sum(regionidx,:)  + Trace_Cortex_Chosen(1,:);

        SPEED = Mice_STRUCT_T(mouseidx).NEURON.Speed(logical(Frame_Inuse_Mouse{mouseidx}))/max(Mice_STRUCT_T(mouseidx).NEURON.Speed(logical(Frame_Inuse_Mouse{mouseidx}))) ;

        hold on;
        plot( SPEED + 1.5*regionidx-1.5, 'LineWidth',1.5,'Color','k' )
        hold on
        text(6900,1.5*regionidx-1,strjoin([string(C_sort(regionidx)) '=' num2str(Neuron_Cortex_Chosen_Num) 'neurons']));


        for trialidx=1:11

            temp = corr( Trace_Cortex_Chosen(result(mouseidx,:)==trialidx)' , SPEED(result(mouseidx,:)==trialidx)' , 'Rows', 'complete' ); % 使用'complete'选项

            R(mouseidx,regionidx,trialidx) = temp;

        end

        hold on;
        plot(45*(Mice_STRUCT_T(mouseidx).NEURON.SHOCK(logical(Frame_Inuse_Mouse{mouseidx}))>0),'k')
        hold on;
        plot(45*(Mice_STRUCT_T(mouseidx).NEURON.CUE(logical(Frame_Inuse_Mouse{mouseidx}))>0),'k')
        hold on;


        plot(60* diff(Trace_unify_SUM_Mice /Trace_unify_NUM_Mice) + 1.5*regionidx-1.51,'k')
        hold on;
    end

    mouseidx
    %    pause

end





f=figure(888)
f.Position = [100 100 200 1000];

for mouseidx=1:5


    for regionidx =1:size(bi_mouse_sorted,2)


        R_mouse_region = squeeze(R(mouseidx,regionidx,:));

        plot(R_mouse_region+regionidx-1, 'LineWidth',1.5,'Color', color_inuse(regionidx,:));
        hold on

        text(0,1*regionidx-1,C_sort(regionidx));
        hold on
        grid on

    end

    % pause
end



R_mouse_mean = squeeze(  mean(R,1)  );

R_mouse_variances = squeeze(  var(R, 0, 1)  );

for regionidx =1:size(bi_mouse_sorted,2)


    plot(R_mouse_mean(regionidx,:)+regionidx-1,'k')
    hold on;
    errorbar(1:11, R_mouse_mean(regionidx,:)+regionidx-1, sqrt(R_mouse_variances(regionidx,:)), 'o');


end






figure(401)

for regionidx =1:size(bi_mouse_sorted,2)

    plot(10* Trace_Cortex_Chosen_Sum(regionidx,:) /Neuron_Cortex_Chosen_Num_Sum(regionidx) + 1*regionidx-1,'Color',color_inuse(regionidx,:));

    hold on
    text(7000,1*regionidx-0.5,C_sort(regionidx));
    hold on

    plot(10* Trace_unify_SUM_Mice /Trace_unify_NUM_Mice + 1*regionidx-1,'k')
    hold on;

end
hold on;
plot(30*(Mice_STRUCT_T(1).NEURON.SHOCK(logical(Frame_Inuse_Mouse{1}))>0),'k')
hold on;
plot(30*(Mice_STRUCT_T(1).NEURON.CUE(logical(Frame_Inuse_Mouse{1}))>0),'k')
hold on;


figure(402)

for regionidx =1:size(bi_mouse_sorted,2)

    plot(10* Trace_Cortex_Chosen_Sum(regionidx,:) /Neuron_Cortex_Chosen_Num_Sum(regionidx) + 1*regionidx-1,'Color',color_inuse(regionidx,:));

    hold on
    text(7000,1*regionidx-0.5,C_sort(regionidx));
    hold on

    plot(10* Trace_unify_SUM_Mice /Trace_unify_NUM_Mice + 1*regionidx-1,'k')
    hold on;

end
hold on;
plot(30*(Mice_STRUCT_T(1).NEURON.SHOCK(logical(Frame_Inuse_Mouse{1}))>0),'k')
hold on;
plot(30*(Mice_STRUCT_T(1).NEURON.CUE(logical(Frame_Inuse_Mouse{1}))>0),'k')
hold on;




%      for i  = 1:numel(miniscopeshockbeginFrame)
%         % 定义时间间隔
%         t1 = miniscopeshockbeginFrame(i)-20; % 开始时间
%         t2 = miniscopeshockendFrame(  i)-20; % 结束时间
%         % 添加阴影
%         px = [t1, t2, t2, t1]; % 阴影的 x 坐标
%         py = [-2, -2, ylimit+2, ylimit+2]; % 阴影的 y 坐标
%         patch(px, py, 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none'); % 添加阴影
%         hold on;
%         % 定义时间间隔
%         t1 = miniscopecuebeginFrame(i); % 开始时间
%         t2 = miniscopecueendFrame(  i); % 结束时间
%         % 添加阴影
%         px = [t1, t2, t2, t1]; % 阴影的 x 坐标
%         py = [-2, -2, ylimit+2, ylimit+2]; % 阴影的 y 坐标
%         patch(px, py, 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none'); % 添加阴影
% 
% 
%     end   