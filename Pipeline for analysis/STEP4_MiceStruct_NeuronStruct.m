%%=========================================================================
%%  ---- This code is for analyzing the fear conditionning experiment 
%%
%%  ---- by Changliang Guo 
%%
%%  ---- 05/15/2024
%%=========================================================================

clc;

clear all;

close all;

%%=========================================================================

mouseidx = 1


for mouseidx =1:5


clearvars -except mouseidx  


Dir                               = 'D:\XM\';

MiceName_T                        = {'FC_M94' ;'FC_M367'; 'FC_M369'; 'FC_M370' ;'FC_M371'};

    MiceName                      = MiceName_T{mouseidx};
    disp(['Processing Mouse: ' MiceName])

path                              = [Dir MiceName '\'];



load(  [path 'shock_cue_frame.mat    ']  );
load(  [path,'mapped_results.mat     ']  );
load(  [path ,'My_WebCam\Behavior.mat']  );
load(  [path 'ROIs_SingleCellComparison_Fig4    ']  );

RAW_ACT                 = RAW_ACT( COMET_Correction.A_good_idx+1 , : ) ;
[RawAct_unify,PS2]      = mapminmax( RAW_ACT,0,1);

RawAct_unify_mean       = mean( RawAct_unify);
RawAct_unify_mean_norm  = mapminmax(RawAct_unify_mean,0,1);

for brainid   = 1:size(COMET_Correction.TraceSUM  ,1)
   
  COMET_Correction.RawAct_uniTraceAVG(brainid,:) = mean( RawAct_unify(find(COMET_Correction.CellRegion_id == brainid),:)  );

end
COMET_Correction.RawAct_unify=RawAct_unify;
COMET_Correction.RawAct_unify_mean=RawAct_unify_mean;
COMET_Correction.RawAct_unify_mean_norm=RawAct_unify_mean_norm;

miniscope_timestamps             = readtable([path 'timeStamps.csv']);
timestamps                       = miniscope_timestamps.TimeStamp_ms_;
timestamps_                      = timestamps(31:end-30);
 


T_spk_diff                       = [0;diff(timestamps_)];
T_spk_diff_Avg                   = mean( T_spk_diff(T_spk_diff>0))/1000;



speed_COMET_smooth               = abs(smoothts(BEHAVIOR.speed_COMET','b',ceil(1/T_spk_diff_Avg)))';
speed_COMET_smooth_7006          = speed_COMET_smooth(31:end-30);
speed_COMET_smooth_7006_binary   = ~isnan(speed_COMET_smooth_7006);
speed_COMET_smooth_7006_binary   = speed_COMET_smooth_7006_binary & (speed_COMET_smooth_7006>=mean(speed_COMET_smooth_7006(speed_COMET_smooth_7006>0)) / 2);
speed_COMET_smooth_7006_interp   = interp1( find(~isnan(speed_COMET_smooth_7006)) , speed_COMET_smooth_7006(~isnan(speed_COMET_smooth_7006)) , 1:1:length(speed_COMET_smooth_7006));   


figure(1001)
plot(speed_COMET_smooth_7006_interp)


%%-------------------------------------------------------------------------


COMET_Correction.CellRegion_id   = atlas_new.valid_id;
COMET_Correction.CellRegion_name = atlas_new.xlabel_annoy;


%%=========================================================================

%%-----------------brain regions atlas-------------------------------------
%%-----------------Correct the atlas mask according to the images online---

figure(100)
    imshow(   atlas_new.pruned_atlas_mask_crop2COMET   )

figure(101)
    imshow(   uint8(label2rgb(atlas_new.clean_top_projection_LR2COMET))   ) 


instrument_size                 = [COMET_Correction.ML_pos_pixel_Scale, COMET_Correction.ML_pos_pixel_Scale];;
instrument_pixel_size           = COMET_Correction.instrument_pixel_size;
CCF_pixel_size                  = atlas_new.CCF_pixel_size;
FOV_size_AP                     = (instrument_size(1) * instrument_pixel_size); %FOV AP size (unit:mm)
FOV_size_ML                     = (instrument_size(2) * instrument_pixel_size); %FOV ML size (unit:mm)
AP_pos_pixel                    = FOV_size_AP / CCF_pixel_size;  %FOV AP size (unit:pixel)
ML_pos_pixel                    = FOV_size_ML / CCF_pixel_size;  %FOV ML size (unit:pixel)
top_projection3                 = atlas_new.clean_top_projection_LR;
P_pos_pixel                     = (instrument_size * instrument_pixel_size) / CCF_pixel_size; 
scale                           = [instrument_size(1)/AP_pos_pixel,instrument_size(2)/ML_pos_pixel];
top_projection3_Scale           = imresize(top_projection3,[scale(1)*size(top_projection3,1), scale(2)*size(top_projection3,2)],"nearest");



figure(102)
    imshow(  label2rgb(top_projection3_Scale),'Colormap',parula  );


labeledImage                    = atlas_new.clean_top_projection_LR2COMET;
angleRad                        = deg2rad(-COMET_Correction.Angle2rot_rightclock);
COMET_cranial_posterior_center  = COMET_Correction.COMET_cranial_posterior_center;
Size_image_Row                  = size(labeledImage,1);    
labeledImage_c                  = imtranslate(double(labeledImage),[0, (COMET_cranial_posterior_center(2)-Size_image_Row/2)],'FillValues',0);
angleRad                        = deg2rad(COMET_Correction.Angle2rot_rightclock);
R                               = [cos(angleRad), -sin(angleRad), 0; sin(angleRad), cos(angleRad), 0; 0, 0, 1];
T1                              = [1, 0, -COMET_cranial_posterior_center(1); 0, 1, -COMET_cranial_posterior_center(2); 0, 0, 1];% 将旋转中心移动到原点
T2                              = [1, 0, COMET_cranial_posterior_center(1); 0, 1, COMET_cranial_posterior_center(2); 0, 0, 1];% 将旋转中心移回原来位置
T                               = T2 * R * T1;% 创建仿射变换矩阵




[height, width, ~]              = size(labeledImage);% 获取图像尺寸
[X, Y]                          = meshgrid(1:width, 1:height);% 创建网格
coords                          = [X(:) Y(:) ones(numel(X), 1)]';
new_coords                      = T * coords;% 应用仿射变换
new_x                           = round(reshape(new_coords(1,:), height, width));% 对新坐标四舍五入并转换为整数
new_y                           = round(reshape(new_coords(2,:), height, width));

labeledImage_c                  = interp2(X, Y, double(labeledImage_c), new_x, new_y,'nearest');
labeledImage_c(isnan(labeledImage_c)) = 0;



addpath([path 'othercolor'],'begin')
Mycolor                           = othercolor('Cat_12');
Mycolor                           = Mycolor(1:1/88*256:256,:);
%%---------------------show the brain region and ID in the figure----------
centroids                   = zeros(length(atlas_new.idLR) , 2);

for grayLevel               = 1:length(atlas_new.idLR) 

    binaryImage             = atlas_new.clean_top_projection_LR2COMET == grayLevel;  

    if sum(binaryImage,"all")

        stats               = regionprops(binaryImage, 'Centroid','Area');    
        temp                = [];

        for idx=1:size(stats,1)
            temp(idx)       = stats(idx).Area;
        end

        idx_max             = find(temp==max(temp));

        if(idx_max)
            centroids(grayLevel, :) = stats(idx_max).Centroid;
        end
       grayLevel
       
    end
 
end
%%---------------------Rotate the figure above clockwirsely by 90 degrees--
centroids_rot90              = centroids;
for idx=1:size(centroids,1)

    centroids_rot90(idx,1)   = 1944-centroids(idx,2);
    centroids_rot90(idx,2)   = centroids(idx,1);

end

centroids_rot90(:,3)         = 30;

atlas_COMET_mask             = (atlas_new.pruned_atlas_mask_crop2COMET>=1)*255;
atlas_COMET_mask_90          = rot90(atlas_COMET_mask,-1);
%%------------------START---Smooth the bourndary
% 形态学操作 - 开运算和闭运算
se                           = strel('disk', 2); % 结构元素，半径为 2 的圆盘
BW_opened                    = imopen(atlas_COMET_mask_90, se); % 开运算
BW_closed                    = imclose(BW_opened, se); % 闭运算
% 骨架化
BW_skeleton                  = bwmorph(BW_closed, 'skel', Inf);
% 高斯滤波器
BW_double                    = double(BW_skeleton); % 转换为 double 类型
BW_filtered                  = imgaussfilt(BW_double, 1); % 高斯滤波，sigma 为 1
% 二值化处理后的图像
atlas_COMET_mask_90_smoothed = imbinarize(BW_filtered);
% 定义结构元素
se                           = strel('disk', 5); % 'disk' 形状的结构元素，半径为 5
% 对图像进行膨胀
atlas_COMET_mask_90_dilate   = imdilate(atlas_COMET_mask_90_smoothed, se);
%%------------------END----Smooth the bourndary

atlas_COMET_mask_90_dilate         = 1 - atlas_COMET_mask_90_dilate;

atlas_COMET_mask_90_dilate_clean   = atlas_COMET_mask_90_dilate;


for idx=1:size(centroids,1)
    if(centroids(idx,1))

        atlas_COMET_mask_90_dilate = insertShape((atlas_COMET_mask_90_dilate),'FilledCircle', round(centroids_rot90(idx,:)), 'Color', Mycolor(idx,:) );
    
    end
end


for idx=1:size(centroids,1)
    if(centroids(idx,1))
        atlas_COMET_mask_90_dilate = insertText(atlas_COMET_mask_90_dilate,round(centroids_rot90(idx,1:2)),num2str(idx),FontSize=30,AnchorPoint='Center',BoxColor='white', BoxOpacity=0,TextColor="black");
    end
end


h=figure(506)
imshow(uint8(atlas_COMET_mask_90_dilate*255));
axis equal
hold on
colormap(Mycolor)
label_names      = string(atlas_new.nameLR);
for idx          = 1:size(label_names,2)
label_names(idx) = strjoin([num2str(idx) '-' label_names(idx)]);
end
colorbar('Ticks',0:1/size(label_names,2):1-1/size(label_names,2),'TickLabels', label_names,'FontSize', 8)
axis tight
axis off

Altas_Image_Regions = frame2im(getframe(h,[0 0 800 800]));
figure(507)
imshow(Altas_Image_Regions);

%%=========================================================================
%% Define interest frames and the reference frames


C_trace                              = COMET_Correction.C_trace;
A_neuron_good_idx                    = COMET_Correction.A_good_idx;

[C_trace_unify_good,PS2]             = mapminmax( C_trace(A_neuron_good_idx+1,:),0,1 );

C_trace_good                         = C_trace(A_neuron_good_idx+1,:);
COMET_Correction.C_trace_unify_good  = C_trace_unify_good;
COMET_Correction.C_trace_good        = C_trace_good;
C_trace_unify_good_diff              = cat( 2 , zeros(size(C_trace_unify_good,1),1) , diff(C_trace_unify_good,1,2) );


%------------------------------------------------------------------------

addpath(genpath([path '\Matlab codes\XM_code\deconv\code\OASIS_matlab']), '-begin');

d1                                   = designfilt('lowpassfir','PassbandFrequency',0.3,'StopbandFrequency',0.4 );

parfor cellno                        = 1:size(C_trace_unify_good,1)   

    C_trace_unify_good_flt(cellno,:) = filtfilt(d1,double(C_trace_unify_good(cellno,:)));
    cellno
end


dt_m                                 = 1;
yshift                               = 1;
idrange                              = 1:101;

figure(103);

for cellno    = idrange


    C_raw_fn  = C_trace_unify_good_flt(cellno,:);
    C_rawn    = C_trace_unify_good(cellno,:);

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

%%-------------------------------------------------------------------------

%------------------------------------------------------------------------
%Adjust the moment of shock: minus by 20 frames


miniscopeshockbeginFrame = miniscopeshockbeginFrame ;
miniscopeshockendFrame   = miniscopeshockendFrame   ;


%%=========================================================================

%%- the frames for spontenous activity

COMET_Correction.sponframe           = zeros(1,size(C_trace_good,2));

spontrial                            = 10;

window_size                          = 200;

for sponid                           = 1:spontrial

    COMET_Correction.sponframe(1,miniscopecuebeginFrame(1)-1-sponid*window_size:miniscopecuebeginFrame(1)-1-(sponid-1)*window_size-1) = spontrial + 1 - sponid;
   
    NumFramespon(sponid)             = sum(COMET_Correction.sponframe==sponid);

    miniscopecuebeginFrame(1)-1-sponid*window_size

    miniscopecuebeginFrame(1)-1-(sponid-1)*window_size-1

end

figure(200)

plot(COMET_Correction.sponframe)

COMET_Correction.sponframe            = zeros(1,size(C_trace_good,2));
COMET_Correction.sponframe(1,50:miniscopecuebeginFrame(1)-1) = 1;
%%- the frames for spontenous activity

%%=========================================================================

%%- the frames for cue activity


miniscopecueendFrame_old               = miniscopecueendFrame;

miniscopecueFrame_max                  = max(miniscopecueendFrame_old - miniscopecuebeginFrame)

miniscopecueendFrame                   = miniscopecuebeginFrame + miniscopecueFrame_max %% make sure each trial has the same number of frames

COMET_Correction.cueframe              = zeros(1,size(C_trace_good,2));

cueid                                  = [];

for i                                  = 1:numel(miniscopecuebeginFrame)
    cueid                              = miniscopecuebeginFrame(i):miniscopecueendFrame(i);
    COMET_Correction.cueframe(1,cueid) = i;
end

figure(201)
plot(COMET_Correction.cueframe)


COMET_Correction.cueframe_2s           = zeros(1,size(C_trace_good,2));
for i                                  = 1:numel(miniscopecuebeginFrame)
    cueid                              = miniscopecuebeginFrame(i):miniscopecuebeginFrame(i)+10;
    COMET_Correction.cueframe_2s(1,cueid) = i;
end

figure(201)
plot(COMET_Correction.cueframe_2s)


%%- the frames for cue activity

%%=========================================================================

%%- the frames for shock activity

COMET_Correction.shockframe                = zeros(1,size(C_trace_good,2));

shockid                                    = [];

for i                                      = 1:numel(miniscopeshockbeginFrame)
    shockid                                = [miniscopeshockbeginFrame(i):miniscopeshockendFrame(i)];
    COMET_Correction.shockframe(1,shockid) = i;
end


figure(202)

    plot(COMET_Correction.shockframe)

%%- the frames for shock activity

%%=========================================================================

%%----- the frames between shock and cues

COMET_Correction.shock2cueframe   = zeros(1,size(C_trace_good,2));

shock2cueid                       = [];

for i                             = 1:numel(miniscopeshockbeginFrame)-1

    shock2cueid                   = [miniscopeshockendFrame(i)+1:miniscopecuebeginFrame(i+1)-1];

    COMET_Correction.shock2cueframe(1,shock2cueid) = i;

end

figure(203)
    plot(COMET_Correction.shock2cueframe)

%%----- the frames between shock and cues

%%=========================================================================

%%----- the frames between shock and cues

COMET_Correction.cue2shockframe  = zeros(1,size(C_trace_good,2));

cue2shockid                      = [];

for i                            = 1:numel(miniscopeshockbeginFrame)
    cue2shockid                  = [miniscopecueendFrame(i)+1:miniscopeshockbeginFrame(i)-1];
    COMET_Correction.cue2shockframe(1,cue2shockid) = i;
end

figure(204)
    plot(COMET_Correction.cue2shockframe)

%%----- the frames between shock and cues

%%=========================================================================

%%----- the frames after the last shock

COMET_Correction.aftershockframe = zeros(1,size(C_trace_good,2));

aftershockid                     = [];

aftershocktrial                  = 3;

window_size                      = 200;

for i                            = 1:aftershocktrial

    aftershockid                 = [miniscopeshockendFrame(end)+window_size*(i-1):miniscopeshockendFrame(end)+window_size*(i)];
    COMET_Correction.aftershockframe(1,aftershockid) = i;

end




  COMET_Correction.aftershockframe = zeros(1,size(C_trace_good,2));
 COMET_Correction.aftershockframe(1,miniscopeshockendFrame(5)+1:miniscopeshockendFrame(5)+400) = 1;
figure(205)
    plot(COMET_Correction.aftershockframe)
%%----- the frames after the last shock

%%=========================================================================

%%----- the frames without cue and shock

COMET_Correction.Other           = ~(COMET_Correction.cueframe + COMET_Correction.shockframe );

figure(206)
    plot(COMET_Correction.Other,'*')

%%----- the frames without cue and shock

%%=========================================================================

%%----- the frames for Cue and others without shock

COMET_Correction.CueplusOther    = COMET_Correction.cueframe*2 + double(COMET_Correction.Other);

figure(207)
plot(COMET_Correction.CueplusOther)


COMET_Correction.Cue2splusOther    = COMET_Correction.cueframe_2s*2 + double(COMET_Correction.Other);

figure(207)
plot(COMET_Correction.Cue2splusOther)
%%----- the frames for Cue and others without shock

%%=========================================================================

%%----- the frames for shock and others without cue

COMET_Correction.ShockplusOther = COMET_Correction.shockframe*2 + double(COMET_Correction.Other);

figure(208)
    plot(COMET_Correction.ShockplusOther)

%%----- the frames for shock and others without cue

%%=========================================================================
COMET_Correction.ShockCueplusOther     = COMET_Correction.shockframe*2 + COMET_Correction.cueframe*2 + double(COMET_Correction.Other);

COMET_Correction.cue2shockplusOther    = COMET_Correction.cue2shockframe*2 + double(COMET_Correction.Other);

COMET_Correction.cue2shockplusOther(COMET_Correction.cue2shockplusOther >2) = COMET_Correction.cue2shockplusOther(COMET_Correction.cue2shockplusOther >2)  -1;

COMET_Correction.cuecue2shockplusOther = COMET_Correction.cueframe*2 + COMET_Correction.cue2shockplusOther;

COMET_Correction.shock2cueplusOther    = COMET_Correction.shock2cueframe*2 + double(COMET_Correction.Other);

COMET_Correction.shock2cueplusOther(COMET_Correction.shock2cueplusOther>1) = COMET_Correction.shock2cueplusOther(COMET_Correction.shock2cueplusOther>1) - 1;

figure(209)
plot( COMET_Correction.sponframe       );
hold on;
plot( COMET_Correction.cueframe        );
hold on;
plot( COMET_Correction.shockframe      );
hold on;
plot( COMET_Correction.shock2cueframe  );
hold on;
plot( COMET_Correction.cue2shockframe  );
hold on;
plot( COMET_Correction.aftershockframe );


figure(2091)
temp1 = COMET_Correction.sponframe;
temp1(temp1>0) = 1;
% plot( temp1 );

temp2=COMET_Correction.cueframe;
temp2(temp2>0) = temp2(temp2>0)*4-2;
temp3=COMET_Correction.cue2shockframe;
temp3(temp3>0) = temp3(temp3>0)*4-1;
temp4=COMET_Correction.shockframe;
temp4(temp4>0) = temp4(temp4>0)*4;
temp5=COMET_Correction.shock2cueframe;
temp5(temp5>0) = temp5(temp5>0)*4+1;
COMET_Correction.trials = temp1 + temp2 + temp3 + temp4+ temp5+COMET_Correction.aftershockframe*21;
plot( COMET_Correction.trials  );



figure(210)
plot( COMET_Correction.ShockCueplusOther );

figure(211)
plot(COMET_Correction.cue2shockplusOther );

figure(212)
plot(COMET_Correction.cuecue2shockplusOther );

figure(213)
plot( COMET_Correction.shock2cueplusOther);

figure(214)
plot(COMET_Correction.shockframe>0 )


temp1                             = COMET_Correction.sponframe';
temp1(temp1>=1)                   = 1;
temp1(temp1<1)                    = 0;
temp1                             = temp1*2;
temp                              = COMET_Correction.shock2cueplusOther';
temp(temp<2)                      = -1;
temp(temp>0)                      = temp(temp>0)+1;

temp                              = temp + temp1;

Behavior_SPON_Trial               = temp + 1;

Behavior_SPON_Trial(Behavior_SPON_Trial>0) = Behavior_SPON_Trial(Behavior_SPON_Trial>0) -1;
Behavior_SPON_Trial = Behavior_SPON_Trial +1;
figure(215)
plot(Behavior_SPON_Trial)

%-----------------------------
COMET_Correction.SponShock2Cue = Behavior_SPON_Trial;
%-----------------------------

temp(temp>0)                      = temp(temp>0)  .*  speed_COMET_smooth_7006_binary(temp>0);
temp                              = temp + 1;
Behavior_SPEED_SponShock2Cue_Trial= temp;

figure(216)
plot(Behavior_SPEED_SponShock2Cue_Trial)

%-----------------------------
COMET_Correction.SPEED_SponShock2Cue = Behavior_SPEED_SponShock2Cue_Trial;
%-----------------------------

temp                              = COMET_Correction.cue2shockplusOther';
temp1 = temp;
temp1(temp1<2) = 0;
temp1(temp1>=2) = 1;
temp2 = temp;
temp2(temp2<2) = 0;

temp1(temp>=2)                    = temp1(temp>=2)  +  speed_COMET_smooth_7006_binary(temp>=2) .* temp2(temp>=2);
temp1(temp1>2) = temp1(temp1>2) - 1;
Behavior_SPEED_Cue2Shock_Trial    = temp1;


figure(217)
plot(Behavior_SPEED_Cue2Shock_Trial)

%-----------------------------
COMET_Correction.SPEED_Cue2Shock = Behavior_SPEED_Cue2Shock_Trial;
%-----------------------------


temp                              = COMET_Correction.CueplusOther';
temp1 = temp;
temp1(temp1<2) = 0;
temp1(temp1>=2) = 1;

temp2 = temp;
temp2(temp2<2) = 0;

base = ones(length(temp));

temp1(temp>=2)              = temp1(temp>=2)  +  speed_COMET_smooth_7006_binary(temp>=2) .* temp2(temp>=2);
temp1(temp1>2)              = temp1(temp1>2) - 1;
Behavior_SPEED_Cue_Trial    = temp1;

figure(218)
plot(Behavior_SPEED_Cue_Trial)
%-----------------------------
COMET_Correction.SPEED_Cue  = Behavior_SPEED_Cue_Trial;
%-----------------------------
%--------------------------------------------------------------------------




Neurons_STRUCT              = {};
Neurons_STRUCT.Trace        = COMET_Correction.C_trace_good;
Neurons_STRUCT.Trace_unify  = COMET_Correction.C_trace_unify_good;
Neurons_STRUCT.Trace_unify_Region_Avg  = COMET_Correction.uniTraceAVG;
Neurons_STRUCT.Trace_unify_Avg  = sum(Neurons_STRUCT.Trace_unify,1)/size(Neurons_STRUCT.Trace_unify,1);
Neurons_STRUCT.Speed        = speed_COMET_smooth_7006';
Neurons_STRUCT.Speed_binary = speed_COMET_smooth_7006_binary';
Neurons_STRUCT.SPON         = COMET_Correction.sponframe;
Neurons_STRUCT.CUE          = COMET_Correction.cueframe;

Neurons_STRUCT.CUE2SHOCK    = COMET_Correction.cue2shockframe;
Neurons_STRUCT.SHOCK        = COMET_Correction.shockframe;

Neurons_STRUCT.LEARN        = COMET_Correction.shock2cueframe;
Neurons_STRUCT.TIME         = timestamps_';
Neurons_STRUCT.CUE_Frame(1,:)   = miniscopecuebeginFrame;
Neurons_STRUCT.CUE_Frame(2,:)   = miniscopecueendFrame;
Neurons_STRUCT.SHOCK_Frame(1,:) = miniscopeshockbeginFrame;
Neurons_STRUCT.SHOCK_Frame(2,:) = miniscopeshockendFrame;
Neurons_STRUCT.CUESHOCK_Frame   = sort( [miniscopecuebeginFrame,miniscopeshockbeginFrame] );
Neurons_STRUCT.A              = COMET_Correction.A_sparse;
Neurons_STRUCT.A_goodidx      = COMET_Correction.A_good_idx;
Neurons_STRUCT.RAW_ACT        = RAW_ACT;
Neurons_STRUCT.RawAct_unify   = COMET_Correction.RawAct_unify;
Neurons_STRUCT.RawAct_uniTraceAVG = COMET_Correction.RawAct_uniTraceAVG;

Neurons_STRUCT.RawAct_unify_mean=COMET_Correction.RawAct_unify_mean;
Neurons_STRUCT.RawAct_unify_mean_norm=COMET_Correction.RawAct_unify_mean_norm;

Mouse_STRUCT                  = {};
Mouse_STRUCT.NAME             = MiceName;
Mouse_STRUCT.NEURON           = Neurons_STRUCT;
Mouse_STRUCT.REGION_ID        = COMET_Correction.Num_neuron';
Mouse_STRUCT.REGION_NAME      = atlas_new.nameLR';
Mouse_STRUCT.REGION_ACT       = mapminmax( COMET_Correction.CortexRegion_uniActivity_raw , 0 , 1 );


Mouse_STRUCT.COORDINATES      = COMET_Correction.N_Rcenter;
Mouse_STRUCT.Neuron_REGION_ID    = COMET_Correction.CellRegion_id';
Mouse_STRUCT.Neuron_REGION_Name  = COMET_Correction.CellRegion_name';
Mouse_STRUCT.CUE              = COMET_Correction.Cell_cueup_idx - COMET_Correction.Cell_cuedown_idx;
Mouse_STRUCT.SHOCK            = COMET_Correction.Cell_Shockup_idx - COMET_Correction.Cell_Shockdown_idx;
Mouse_STRUCT.CENTROIDS90      = centroids_rot90;
Mouse_STRUCT.ATLAS90          = atlas_COMET_mask_90_dilate_clean;
Mouse_STRUCT.ATLAS2COMET      = atlas_new.clean_top_projection_LR2COMET;
Mouse_STRUCT.ATLAS2SCALE      = top_projection3_Scale;
Mouse_STRUCT.ATLASFIG         = Altas_Image_Regions;


save([Dir  ,'/' MiceName '_Mouse_STRUCT.mat'], 'Mouse_STRUCT', '-v7.3');

%%=========================================================================

% pause 

end


