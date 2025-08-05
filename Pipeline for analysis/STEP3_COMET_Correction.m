%%=========================================================================
%%
%%  ----This code is for analyzing the fear conditionning experiment 
%%
%%  ----by Changliang Guo 
%%
%%  ----05/15/2024
%%
%%=========================================================================

clc;

clear;

close all;

%%=========================================================================

path                              = 'D:\XM\FC_M94\';

% load([path 'AMF_despeckle_MC_denoised_8bit_caiman_result_handcured.mat']);

load(  [path 'shock_cue_frame.mat    '] );

load(  [path,'mapped_results.mat     '] );

load(  [path ,'My_WebCam\Behavior.mat'] );

miniscope_timestamps             = readtable([path 'timeStamps.csv']);

timestamps                       = miniscope_timestamps.TimeStamp_ms_;

timestamps_                      = timestamps(31:end-30);

T_spk_diff                       = [0;diff(timestamps_)];
%%-------------------------------------------------------------------------

T_spk_diff_Avg                   = mean( T_spk_diff(T_spk_diff>0))/1000;

speed_COMET_smooth               = abs(smoothts(BEHAVIOR.speed_COMET','b',ceil(1/T_spk_diff_Avg)))';

speed_COMET_smooth_7006          = speed_COMET_smooth(31:end-30);

speed_COMET_smooth_7006_binary   = ~isnan(speed_COMET_smooth_7006);

speed_COMET_smooth_7006_binary   = speed_COMET_smooth_7006_binary & (speed_COMET_smooth_7006>=mean(speed_COMET_smooth_7006(speed_COMET_smooth_7006>0)) / 1);

figure(1000)

    plot(speed_COMET_smooth_7006)
    hold on;
    plot( speed_COMET_smooth_7006_binary*40   ) 
    hold on;
    bar(speed_COMET_smooth_7006_binary*40 )

%%-------------------------------------------------------------------------

COMET_Correction.CellRegion_id   = atlas_new.valid_id;

COMET_Correction.CellRegion_name = atlas_new.xlabel_annoy;


%%=========================================================================

%%-----------------brain regions atlas-------------------------------------
%%-----------------Correct the atlas mask according to the images online---


figure(100)
    imshow(atlas_new.pruned_atlas_mask_crop2COMET)

figure(101)
    imshow(uint8(label2rgb(atlas_new.clean_top_projection_LR2COMET)))


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
    imshow( label2rgb(top_projection3_Scale),'Colormap',parula );


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


%%=========================================================================

%% Define interest frames and the reference frames

% spk归一

% [S_trace_unify_good,PS] = mapminmax(s_2d,0,1);

C_trace                              = COMET_Correction.C_trace;

A_neuron_good_idx                    = COMET_Correction.A_good_idx;

[C_trace_unify_good,PS2]             = mapminmax(C_trace(A_neuron_good_idx+1,:),0,1);

C_trace_good                         = C_trace(A_neuron_good_idx+1,:);

COMET_Correction.C_trace_unify_good  = C_trace_unify_good;

COMET_Correction.C_trace_good        = C_trace_good;



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

    shock2cueid                   = [miniscopeshockendFrame(i):miniscopecuebeginFrame(i+1)];

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
    cue2shockid                  = [miniscopecueendFrame(i):miniscopeshockbeginFrame(i)];
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

COMET_Correction.CueplusOther   = COMET_Correction.cueframe*2 + double(COMET_Correction.Other);

figure(207)
    plot(COMET_Correction.CueplusOther)

%%----- the frames for Cue and others without shock

%%=========================================================================

%%----- the frames for shock and others without cue

COMET_Correction.ShockplusOther = COMET_Correction.shockframe*2 + double(COMET_Correction.Other);

figure(208)
    plot(COMET_Correction.ShockplusOther)

%%----- the frames for shock and others without cue

%%=========================================================================


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
    

% for grpidx=1:100
% 
%     figure(210)    
% 
%     
%     for tridx=1:100
%         plot(C_trace_unify_good(tridx+100*(grpidx-1),:)+tridx);
%         xlim([0 7009])
%         ylim([0 101])
%         hold on
%     end
% 
%     for i  = 1:numel(miniscopecuebeginFrame)
%         % 定义时间间隔
%         t1 = miniscopecuebeginFrame(i); % 开始时间
%         t2 = miniscopecueendFrame(i); % 结束时间
%         % 添加阴影
%         px = [t1, t2, t2, t1]; % 阴影的 x 坐标
%         py = [-2, -2, 100+2, 100+2]; % 阴影的 y 坐标
%         patch(px, py, 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none'); % 添加阴影
% 
%      
%     end
%     grpidx
% %     pause
%     hold off
% 
% end


%%=========================================================================
%%=========================================================================

%%===Analyze cue related brain regions

Cortex4Cue                        = COMET_Correction.CortexRegion_uniActivity_raw(:,COMET_Correction.CueplusOther>0);
Frameidx_CueplusOther             = COMET_Correction.CueplusOther(COMET_Correction.CueplusOther~=0);
figure(300)
plot(Frameidx_CueplusOther)


ns                                = size(Cortex4Cue,1);
numShuffles                       = 500;
numShuffles                       = numShuffles+1;

NoFrames                          = size(Cortex4Cue,2);
shift                             = floor( rand(numShuffles,1) *NoFrames);
shift(1,1)                        = 0;
indexp                            = 1:NoFrames;


T_spk_diff_CueOther               = T_spk_diff(COMET_Correction.CueplusOther>0);
T_spk_diff_cue                    = T_spk_diff(COMET_Correction.cueframe>0  );
T_spk_diff_Other                  = T_spk_diff(COMET_Correction.Other       );

T_spk_diff_cue_sum                = sum(T_spk_diff_cue);
T_spk_diff_Other_sum              = sum(T_spk_diff_Other);
FR_other                          = []

parfor idxnumShf                  = 1:numShuffles

    indexs2                       = mod(indexp+shift(idxnumShf), NoFrames);
    indexs2(indexs2==0)           = NoFrames;
    Cortex4Cue_shf                = Cortex4Cue(:,indexs2);

    Cortex4Cue_shf_cue            = Cortex4Cue_shf(:,Frameidx_CueplusOther>=2);
    Cortex4Cue_shf_cue_T1         = Cortex4Cue_shf(:,Frameidx_CueplusOther==2);
    Cortex4Cue_shf_cue_T2         = Cortex4Cue_shf(:,Frameidx_CueplusOther==4);
    Cortex4Cue_shf_cue_T3         = Cortex4Cue_shf(:,Frameidx_CueplusOther==6);
    Cortex4Cue_shf_cue_T4         = Cortex4Cue_shf(:,Frameidx_CueplusOther==8);
    Cortex4Cue_shf_cue_T5         = Cortex4Cue_shf(:,Frameidx_CueplusOther==10);
    Cortex4Cue_shf_other          = Cortex4Cue_shf(:,Frameidx_CueplusOther==1);
    Cortex4Cue_shf_cue_sum        = sum(Cortex4Cue_shf_cue,2);
    Cortex4Cue_shf_cue_T1_sum     = sum(Cortex4Cue_shf_cue_T1,2);
    Cortex4Cue_shf_cue_T2_sum     = sum(Cortex4Cue_shf_cue_T2,2);
    Cortex4Cue_shf_cue_T3_sum     = sum(Cortex4Cue_shf_cue_T3,2);
    Cortex4Cue_shf_cue_T4_sum     = sum(Cortex4Cue_shf_cue_T4,2);
    Cortex4Cue_shf_cue_T5_sum     = sum(Cortex4Cue_shf_cue_T5,2);
    Cortex4Cue_shf_other_sum      = sum(Cortex4Cue_shf_other,2);

    FR_Cortex4Cue(   idxnumShf,:) = Cortex4Cue_shf_cue_sum/T_spk_diff_cue_sum;
    FR_Cortex4Cue_T1(idxnumShf,:) = Cortex4Cue_shf_cue_T1_sum/T_spk_diff_cue_sum * 5;
    FR_Cortex4Cue_T2(idxnumShf,:) = Cortex4Cue_shf_cue_T2_sum/T_spk_diff_cue_sum * 5;
    FR_Cortex4Cue_T3(idxnumShf,:) = Cortex4Cue_shf_cue_T3_sum/T_spk_diff_cue_sum * 5;
    FR_Cortex4Cue_T4(idxnumShf,:) = Cortex4Cue_shf_cue_T4_sum/T_spk_diff_cue_sum * 5;
    FR_Cortex4Cue_T5(idxnumShf,:) = Cortex4Cue_shf_cue_T5_sum/T_spk_diff_cue_sum * 5;
    FR_other(idxnumShf,:)         = Cortex4Cue_shf_other_sum/T_spk_diff_Other_sum;

    idxnumShf

end

FR_Cortex4Cue_other               = FR_Cortex4Cue - FR_other;
FR_Cortex4Cue_other_Real          = FR_Cortex4Cue_other(1,:);
infoP_Cortex4Cue                  = sum(repmat(FR_Cortex4Cue_other_Real,numShuffles-1,1) >= FR_Cortex4Cue_other(2:end,:),1)./(numShuffles-1);

FR_Cortex4Cue_T1_other            = FR_Cortex4Cue_T1 - FR_other;
FR_Cortex4Cue_T1_other_Real       = FR_Cortex4Cue_T1_other(1,:);
infoP_Cortex4Cue_T1               = sum(repmat(FR_Cortex4Cue_T1_other_Real,numShuffles-1,1) >= FR_Cortex4Cue_T1_other(2:end,:),1)./(numShuffles-1);

FR_Cortex4Cue_T2_other            = FR_Cortex4Cue_T2 - FR_other;
FR_Cortex4Cue_T2_other_Real       = FR_Cortex4Cue_T2_other(1,:);
infoP_Cortex4Cue_T2               = sum(repmat(FR_Cortex4Cue_T2_other_Real,numShuffles-1,1) >= FR_Cortex4Cue_T2_other(2:end,:),1)./(numShuffles-1);

FR_Cortex4Cue_T3_other            = FR_Cortex4Cue_T3 - FR_other;
FR_Cortex4Cue_T3_other_Real       = FR_Cortex4Cue_T3_other(1,:);
infoP_Cortex4Cue_T3               = sum(repmat(FR_Cortex4Cue_T3_other_Real,numShuffles-1,1) >= FR_Cortex4Cue_T3_other(2:end,:),1)./(numShuffles-1);

FR_Cortex4Cue_T4_other            = FR_Cortex4Cue_T4 - FR_other;
FR_Cortex4Cue_T4_other_Real       = FR_Cortex4Cue_T4_other(1,:);
infoP_Cortex4Cue_T4               = sum(repmat(FR_Cortex4Cue_T4_other_Real,numShuffles-1,1) >= FR_Cortex4Cue_T4_other(2:end,:),1)./(numShuffles-1);

FR_Cortex4Cue_T5_other            = FR_Cortex4Cue_T5 - FR_other;
FR_Cortex4Cue_T5_other_Real       = FR_Cortex4Cue_T5_other(1,:);
infoP_Cortex4Cue_T5               = sum(repmat(FR_Cortex4Cue_T5_other_Real,numShuffles-1,1) >= FR_Cortex4Cue_T5_other(2:end,:),1)./(numShuffles-1);

infoP_Cortex4Cue_T                = []
infoP_Cortex4Cue_T(1,:)           = infoP_Cortex4Cue;     
infoP_Cortex4Cue_T(2,:)           = infoP_Cortex4Cue_T1;    
infoP_Cortex4Cue_T(3,:)           = infoP_Cortex4Cue_T2;    
infoP_Cortex4Cue_T(4,:)           = infoP_Cortex4Cue_T3;    
infoP_Cortex4Cue_T(5,:)           = infoP_Cortex4Cue_T4;    
infoP_Cortex4Cue_T(6,:)           = infoP_Cortex4Cue_T5;    

Num_Cortex4Cue_095                = sum(infoP_Cortex4Cue_T>=0.95,1)
Num_Cortex4Cue_005                = sum(infoP_Cortex4Cue_T<=0.05,1)

infoP_Cortex4Cue_T                = infoP_Cortex4Cue_T';
%%=========================================================================

addpath([path 'othercolor'],'begin')
Mycolor                           = othercolor('Cat_12');
Mycolor                           = Mycolor(1:1/88*256:256,:);

figure(300)
for idx=1:length(atlas_new.idLR) 
    idx
imshow(label2rgb(atlas_new.clean_top_projection_LR2COMET==idx))

% pause
end

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
%%---------------------show the brain region and ID in the figure----------

figure(301)
    imshow(  label2rgb(atlas_new.clean_top_projection_LR2COMET, Mycolor  ) )
    hold on;
    for idx=1:size(centroids,1)
        if(centroids(idx,1))
            text(centroids(idx,1)+0.15,centroids(idx,2),num2str(idx), 'Color', 'black')
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

figure(503)
imshow(  uint8(atlas_COMET_mask_90/1)  );
hold on;

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

figure(504)
imshow((atlas_COMET_mask_90_dilate/1));

atlas_COMET_mask_90_dilate         = 1 - atlas_COMET_mask_90_dilate;

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

figure(506)
imshow(uint8(atlas_COMET_mask_90_dilate*255));
axis equal
hold on
colormap(Mycolor)
label_names      = string(atlas_new.nameLR);
for idx          = 1:size(label_names,2)
label_names(idx) = strjoin([num2str(idx) '-' label_names(idx)]);
end
colorbar('Ticks',0:1/size(label_names,2):1-1/size(label_names,2),'TickLabels', label_names)
axis tight
axis off

%%-------------------------------------------------------------------------
%%---------------------Rotate the figure above clockwirsely by 90 degrees--

% imwrite(uint8(labeledImage_update_rgb),'brain_atlas_color.png')

%%=========================================================================
%%=========================================================================
Mycolor(89,:) = [1 1 1]
%%---------------------------
figure(600)
plot(  infoP_Cortex4Cue_T(sum(COMET_Correction.CortexRegion_uniActivity_raw,2)>0,1)  )

atlas_COMET_mask_90_dilate         = imdilate(atlas_COMET_mask_90_smoothed, se);
atlas_COMET_mask_90_dilate         = 1 - atlas_COMET_mask_90_dilate;

for idx=1:size(centroids,1)

    idx
    flag = (  centroids_rot90(idx,1) & centroids_rot90(idx,2)  )
    if(flag&&(infoP_Cortex4Cue_T(idx,1)>0.0))

        atlas_COMET_mask_90_dilate = insertShape((atlas_COMET_mask_90_dilate),'FilledCircle', round(centroids_rot90(idx,:)), 'Color', Mycolor(idx,:) );
        atlas_COMET_mask_90_dilate = insertText(atlas_COMET_mask_90_dilate,round(centroids_rot90(idx,1:2)),num2str(idx),FontSize=30,AnchorPoint='Center',BoxColor='white', BoxOpacity=0,TextColor="black");
        figure(601)
        imshow(uint8(atlas_COMET_mask_90_dilate*255));
%         pause

    end

end

%%---------------------------

figure(601)
imshow(uint8(atlas_COMET_mask_90_dilate*255));
hold on
colormap(Mycolor)
label_names      = string(atlas_new.nameLR);
for idx          = 1:size(label_names,2)
label_names(idx) = strjoin([num2str(idx) '-' label_names(idx)]);
end
colorbar('Ticks',0:1/size(label_names,2):1-1/size(label_names,2),'TickLabels', label_names)
axis tight
axis off

%%---------------------------

%%=========================================================================

NumY_show                          = 2000
for Cue_Trial                      = 1:size(infoP_Cortex4Cue_T,2)

    [infoP_Cortex4Cue_T_sort  infoP_Cortex4cue_T_Index(:,Cue_Trial)]  = sortrows(infoP_Cortex4Cue_T,-Cue_Trial);
    
    Cortex_cueup(:,Cue_Trial)      = (infoP_Cortex4Cue_T_sort(:,Cue_Trial)>=0.6);
    Cortex_cueup_Num(Cue_Trial)    = sum(double(Cortex_cueup(:,Cue_Trial) ));
       
    Cortex_cuedown(:,Cue_Trial)    = (infoP_Cortex4Cue_T_sort(:,Cue_Trial)<=0.05);
    Cortex_cuedown_Num(Cue_Trial)  = sum(double(Cortex_cuedown(:,Cue_Trial) ));
        
    Cortex_cuebetween(:,Cue_Trial) = (infoP_Cortex4Cue_T_sort(:,Cue_Trial)<0.6)&(infoP_Cortex4Cue_T_sort(:,Cue_Trial)>0.05);
        
    CortexRegion_uniActivity_raw_good       = COMET_Correction.CortexRegion_uniActivity_raw(infoP_Cortex4cue_T_Index(:,Cue_Trial),:);
    CortexRegion_uniActivity_raw_cueup      = CortexRegion_uniActivity_raw_good( Cortex_cueup(:,Cue_Trial),     :);
    CortexRegion_uniActivity_raw_cuedown    = CortexRegion_uniActivity_raw_good( Cortex_cuedown(:,Cue_Trial),   :);
    CortexRegion_uniActivity_raw_cuebetween = CortexRegion_uniActivity_raw_good( Cortex_cuebetween(:,Cue_Trial),:);
    
    
    figure(1000 + Cue_Trial)
   
    for i = 1:min(NumY_show, size(CortexRegion_uniActivity_raw_cueup,1))
        plot(CortexRegion_uniActivity_raw_cueup(i,:)*1+i,'k');
        hold on;
    end
    axis([1, size(CortexRegion_uniActivity_raw_cueup,2), -2, min(NumY_show, size(CortexRegion_uniActivity_raw_cueup,1))+2]);
    
    for i  = 1:numel(miniscopecuebeginFrame)
        % 定义时间间隔
        t1 = miniscopecuebeginFrame(i); % 开始时间
        t2 = miniscopecueendFrame(i); % 结束时间
        % 添加阴影
        px = [t1, t2, t2, t1]; % 阴影的 x 坐标
        py = [-2, -2, min(NumY_show, size(CortexRegion_uniActivity_raw_cueup,1))+2, min(NumY_show, size(CortexRegion_uniActivity_raw_cueup,1))+2]; % 阴影的 y 坐标
        patch(px, py, 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none'); % 添加阴影
    
        line([t1 t1], [-2 min(NumY_show, size(CortexRegion_uniActivity_raw_cueup,1))+2], 'Color', 'r', 'LineStyle', '--'); % 添加开始时间的竖直线
        line([t2 t2], [-2 min(NumY_show, size(CortexRegion_uniActivity_raw_cueup,1))+2], 'Color', 'r', 'LineStyle', '--'); % 添加结束时间的竖直线
    
    end
    title('The first 1000 cue up cell trace')


    figure(1100 + Cue_Trial)

    for i = 1:min(NumY_show, size(CortexRegion_uniActivity_raw_cuedown,1))
        plot(CortexRegion_uniActivity_raw_cuedown(i,:)*1+i,'k');
        hold on;
    end
    axis([1, size(CortexRegion_uniActivity_raw_cuedown,2), -2, min(NumY_show, size(CortexRegion_uniActivity_raw_cuedown,1))+2]);
    
    for i  = 1:numel(miniscopecuebeginFrame)
        % 定义时间间隔
        t1 = miniscopecuebeginFrame(i); % 开始时间
        t2 = miniscopecueendFrame(i); % 结束时间
        % 添加阴影
        px = [t1, t2, t2, t1]; % 阴影的 x 坐标
        py = [-2, -2, min(NumY_show, size(CortexRegion_uniActivity_raw_cuedown,1))+2, min(NumY_show, size(CortexRegion_uniActivity_raw_cuedown,1))+2]; % 阴影的 y 坐标
        patch(px, py, 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none'); % 添加阴影
    
        line([t1 t1], [-2 min(NumY_show, size(CortexRegion_uniActivity_raw_cuedown,1))+2], 'Color', 'r', 'LineStyle', '--'); % 添加开始时间的竖直线
        line([t2 t2], [-2 min(NumY_show, size(CortexRegion_uniActivity_raw_cuedown,1))+2], 'Color', 'r', 'LineStyle', '--'); % 添加结束时间的竖直线
    end

    title('The first 1000 cue down cell trace')

Cue_Trial
    
end 

%%=========================================================================
%%=========================================================================

%%---------START--Plot the brain region activity---------

CortexRegion_uniActivity_raw      = COMET_Correction.CortexRegion_uniActivity_raw;
CortexRegion_uniActivity_raw_norm = (CortexRegion_uniActivity_raw-min(CortexRegion_uniActivity_raw,[],2))./repmat(max(CortexRegion_uniActivity_raw,[],2)-min(CortexRegion_uniActivity_raw,[],2),1,size(CortexRegion_uniActivity_raw,2));

flag                              = sum(CortexRegion_uniActivity_raw,2);
for brainid                       = 1:size(atlas_new.nameLR,2)    
    
    if( flag(brainid) )
        figure(700)
%         plot(CortexRegion_uniActivity_raw(size(atlas_new.nameLR,2)-brainid+1,:)*1,'color',Mycolor(floor((size(atlas_new.nameLR,2)-brainid+1)),:));
        plot(CortexRegion_uniActivity_raw(brainid,:)*1,'color','black');
        title(['Brain region: ' num2str(brainid) '-Name: ' atlas_new.nameLR{1,brainid}])
%         axis tight
%         axis off
%         hold on;
    end  
% set(gcf, 'color', 'none');   
% set(gca, 'color', 'none'); 
brainid   
% pause
end

figure(701)
for brainid                       = 1:size(atlas_new.nameLR,2)
    subplot(size(atlas_new.nameLR,2), 1, brainid)
    plot(CortexRegion_uniActivity_raw(brainid,:),'color',Mycolor(floor((size(atlas_new.nameLR,2)-brainid+1)),:));
    axis tight
    axis off
    brainid
end

figure(702)
for brainid                       = 1:size(atlas_new.nameLR,2)
    
    flag = sum(centroids_rot90(idx,:))
    if( flag )
        plot(CortexRegion_uniActivity_raw_norm(brainid,:)+brainid,'color',Mycolor(floor((size(atlas_new.nameLR,2)-brainid+1)),:));
        hold on;
    end

end

%%---------END--Plot the brain region activity---------

%%=========================================================================
%%=========================================================================

%%------------- Correlation matrix for cue
figure(800)
plot(COMET_Correction.sponframe)
CortexRegion_uniActivity_raw_norm_spon = []
temp=[]
for idx              = 1:max(unique(COMET_Correction.sponframe))

    for brianidx     = 1: size(atlas_new.nameLR,2)
        
        CortexRegion_uniActivity_raw_norm_spon(brianidx,:,idx) = CortexRegion_uniActivity_raw_norm(brianidx,COMET_Correction.sponframe==idx);

    end

idx
end

for idx              = 1:max(unique(COMET_Correction.sponframe))

    temp             = squeeze( CortexRegion_uniActivity_raw_norm_spon(:,:,idx) );
    Correlation_Spon(:,:,idx) = corrcoef(temp');

end


figure(701)
tiledlayout(2,5, 'Padding', 'none', 'TileSpacing', 'compact');

for column = 1:10
    subplot(2,5,column)
    heatmap(squeeze(Correlation_Spon(:,:,column)), 'Colormap', parula);
    colorbar off
    grid off
end

%%=========================================================================
%%-------------
figure(801)
plot(COMET_Correction.cueframe)

CortexRegion_uniActivity_raw_norm_cue = []
%%------------------------------------------the number of each cue frame is
%%slightly different, to make the matrix inconsistant. 

temp                 = []

for idx              = 1:max(unique(COMET_Correction.cueframe))

    for brianidx     = 1: size(atlas_new.nameLR,2)
        
%         CortexRegion_uniActivity_raw_norm_cue(brianidx,:,idx) = CortexRegion_uniActivity_raw_norm(brianidx,COMET_Correction.cueframe==idx);


        CortexRegion_uniActivity_raw_norm_cue{idx}(brianidx,:) = CortexRegion_uniActivity_raw_norm(brianidx,COMET_Correction.cueframe==idx);

    end

idx
end

%%------------------------------------------the number of each cue frame is
%%slightly different, to make the matrix inconsistant. 

for idx              = 1:max(unique(COMET_Correction.cueframe))

    temp             = squeeze( CortexRegion_uniActivity_raw_norm_cue{idx} );
    Correlation_cue(:,:,idx) = corrcoef(temp');

end

figure(801)
tiledlayout(2,5, 'Padding', 'none', 'TileSpacing', 'compact');

for column=1:max(unique(COMET_Correction.cueframe))
    subplot(1,5,column)
    heatmap(squeeze(Correlation_cue(:,:,column)), 'Colormap', parula);
    colorbar off
    grid off
end


temp                 = []

CortexRegion_uniActivity_raw_norm_noNAN = CortexRegion_uniActivity_raw_norm( all(~isnan(CortexRegion_uniActivity_raw_norm),2), : );
CortexRegion_uniActivity_raw_norm_noNAN_cue=[];
for idx              = 1:max(unique(COMET_Correction.cueframe));

    for brianidx     = 1: size(CortexRegion_uniActivity_raw_norm_noNAN,1)
        
%         CortexRegion_uniActivity_raw_norm_noNAN_cue(brianidx,:,idx) = CortexRegion_uniActivity_raw_norm_noNAN(brianidx,COMET_Correction.cueframe==idx);

        CortexRegion_uniActivity_raw_norm_noNAN_cue{idx}(brianidx,:) = CortexRegion_uniActivity_raw_norm_noNAN(brianidx,COMET_Correction.cueframe==idx);
    end

idx
end

Correlation_cue_noNAN = []
for idx               = 1:max(unique(COMET_Correction.cueframe))

    temp = squeeze( CortexRegion_uniActivity_raw_norm_noNAN_cue{idx} );
    Correlation_cue_noNAN(:,:,idx) = corrcoef(temp');

end


figure(802)
tiledlayout(1,5, 'Padding', 'none', 'TileSpacing', 'compact');

for column=1:max(unique(COMET_Correction.cueframe))
    subplot(1,5,column)
    heatmap(squeeze(Correlation_cue_noNAN(:,:,column)), 'Colormap', parula);
    colorbar off
    grid off
end


%%=========================================================================
%%=========================================================================
%%------------- Correlation matrix for shock

figure(900)
plot( COMET_Correction.shockframe )

CortexRegion_uniActivity_raw_norm_shock = []

temp                 = []
for idx              = 1:max(unique(COMET_Correction.shockframe))

    for brianidx     = 1: size(atlas_new.nameLR,2)
        
        CortexRegion_uniActivity_raw_norm_shock{idx}(brianidx,:) = CortexRegion_uniActivity_raw_norm(brianidx,COMET_Correction.shockframe==idx);

    end

idx
end

for idx              = 1:max(unique(COMET_Correction.shockframe))

    temp             = squeeze( CortexRegion_uniActivity_raw_norm_shock{idx} );
    Correlation_shock(:,:,idx) = corrcoef(temp');

end


figure(901)
tiledlayout(2,5, 'Padding', 'none', 'TileSpacing', 'compact');

for column=1:max(unique(COMET_Correction.shockframe))
    subplot(1,5,column)
    heatmap(squeeze(Correlation_shock(:,:,column)), 'Colormap', parula);
    colorbar off
    grid off
end


temp                 = []

CortexRegion_uniActivity_raw_norm_noNAN       = CortexRegion_uniActivity_raw_norm( all(~isnan(CortexRegion_uniActivity_raw_norm),2), : );
CortexRegion_uniActivity_raw_norm_noNAN_shock = [];
for idx               = 1:max(unique(COMET_Correction.shockframe));

    for brianidx      = 1: size(CortexRegion_uniActivity_raw_norm_noNAN,1)
         
        CortexRegion_uniActivity_raw_norm_noNAN_shock{idx}(brianidx,:) = CortexRegion_uniActivity_raw_norm_noNAN(brianidx,COMET_Correction.shockframe==idx);

    end

idx
end

Correlation_shock_noNAN = []
for idx               = 1:max(unique(COMET_Correction.shockframe))

    temp = squeeze( CortexRegion_uniActivity_raw_norm_noNAN_shock{idx} );
    Correlation_shock_noNAN(:,:,idx) = corrcoef(temp');

end

figure(902)
tiledlayout(1,5, 'Padding', 'none', 'TileSpacing', 'compact');



for column           = 1:max(unique(COMET_Correction.shockframe))
    subplot(1,5,column)
    heatmap(squeeze(Correlation_shock_noNAN(:,:,column)), 'Colormap', parula);
    colorbar off
    grid off
end


%%=========================================================================

%%===Analyze cue related neurons


    Trace4Cue                     = C_trace_unify_good(:,COMET_Correction.CueplusOther>0);
    Frameidx_CueplusOther         = COMET_Correction.CueplusOther(COMET_Correction.CueplusOther~=0);
   
    figure(600)
    plot(Frameidx_CueplusOther)
    
    ns                            = size(Trace4Cue,1);
    numShuffles                   = 500;
    numShuffles                   = numShuffles+1;
    
    NoFrames                      = size(Trace4Cue,2);
    shift                         = floor( rand(numShuffles,1) *NoFrames);
    shift(1,1)                    = 0;
    indexp                        = 1:NoFrames;
    
    T_spk_diff                    = [0;diff(timestamps_)];
    
    T_spk_diff_CueOther           = T_spk_diff(COMET_Correction.CueplusOther>0);
    T_spk_diff_cue                = T_spk_diff(COMET_Correction.cueframe>0  );
    T_spk_diff_Other              = T_spk_diff(COMET_Correction.Other       );
    
    T_spk_diff_cue_sum            = sum(T_spk_diff_cue);
    T_spk_diff_Other_sum          = sum(T_spk_diff_Other);
    FR_other=[]

    parfor idxnumShf              = 1:numShuffles
    
        indexs2                   = mod(indexp+shift(idxnumShf), NoFrames);
        indexs2(indexs2==0)       = NoFrames;
        Trace4Cue_shf             = Trace4Cue(:,indexs2);
    
        Trace4Cue_shf_cue         = Trace4Cue_shf(:,Frameidx_CueplusOther>=2);
        Trace4Cue_shf_cue_T1      = Trace4Cue_shf(:,Frameidx_CueplusOther==2);
        Trace4Cue_shf_cue_T2      = Trace4Cue_shf(:,Frameidx_CueplusOther==4);
        Trace4Cue_shf_cue_T3      = Trace4Cue_shf(:,Frameidx_CueplusOther==6);
        Trace4Cue_shf_cue_T4      = Trace4Cue_shf(:,Frameidx_CueplusOther==8);
        Trace4Cue_shf_cue_T5      = Trace4Cue_shf(:,Frameidx_CueplusOther==10);
        Trace4Cue_shf_other       = Trace4Cue_shf(:,Frameidx_CueplusOther==1);
        Trace4Cue_shf_cue_sum     = sum(Trace4Cue_shf_cue,2);
        Trace4Cue_shf_cue_T1_sum  = sum(Trace4Cue_shf_cue_T1,2);
        Trace4Cue_shf_cue_T2_sum  = sum(Trace4Cue_shf_cue_T2,2);
        Trace4Cue_shf_cue_T3_sum  = sum(Trace4Cue_shf_cue_T3,2);
        Trace4Cue_shf_cue_T4_sum  = sum(Trace4Cue_shf_cue_T4,2);
        Trace4Cue_shf_cue_T5_sum  = sum(Trace4Cue_shf_cue_T5,2);
        Trace4Cue_shf_other_sum   = sum(Trace4Cue_shf_other,2);
    
        FR_cue(idxnumShf,:)       = Trace4Cue_shf_cue_sum/T_spk_diff_cue_sum;
        FR_cue_T1(idxnumShf,:)    = Trace4Cue_shf_cue_T1_sum/T_spk_diff_cue_sum * 5;
        FR_cue_T2(idxnumShf,:)    = Trace4Cue_shf_cue_T2_sum/T_spk_diff_cue_sum * 5;
        FR_cue_T3(idxnumShf,:)    = Trace4Cue_shf_cue_T3_sum/T_spk_diff_cue_sum * 5;
        FR_cue_T4(idxnumShf,:)    = Trace4Cue_shf_cue_T4_sum/T_spk_diff_cue_sum * 5;
        FR_cue_T5(idxnumShf,:)    = Trace4Cue_shf_cue_T5_sum/T_spk_diff_cue_sum * 5;
        FR_other(idxnumShf,:)     = Trace4Cue_shf_other_sum/T_spk_diff_Other_sum;
    
        idxnumShf
    
    end
    
    FR_cue_other                  = FR_cue - FR_other;
    FR_cue_other_Real             = FR_cue_other(1,:);
    infoP_cue                     = sum(repmat(FR_cue_other_Real,numShuffles-1,1) >= FR_cue_other(2:end,:),1)./(numShuffles-1);
    
    FR_cue_T1_other               = FR_cue_T1 - FR_other;
    FR_cue_T1_other_Real          = FR_cue_T1_other(1,:);
    infoP_cue_T1                  = sum(repmat(FR_cue_T1_other_Real,numShuffles-1,1) >= FR_cue_T1_other(2:end,:),1)./(numShuffles-1);
    
    FR_cue_T2_other               = FR_cue_T2 - FR_other;
    FR_cue_T2_other_Real          = FR_cue_T2_other(1,:);
    infoP_cue_T2                  = sum(repmat(FR_cue_T2_other_Real,numShuffles-1,1) >= FR_cue_T2_other(2:end,:),1)./(numShuffles-1);
    
    FR_cue_T3_other               = FR_cue_T3 - FR_other;
    FR_cue_T3_other_Real          = FR_cue_T3_other(1,:);
    infoP_cue_T3                  = sum(repmat(FR_cue_T3_other_Real,numShuffles-1,1) >= FR_cue_T3_other(2:end,:),1)./(numShuffles-1);
    
    FR_cue_T4_other               = FR_cue_T4 - FR_other;
    FR_cue_T4_other_Real          = FR_cue_T4_other(1,:);
    infoP_cue_T4                  = sum(repmat(FR_cue_T4_other_Real,numShuffles-1,1) >= FR_cue_T4_other(2:end,:),1)./(numShuffles-1);
    
    FR_cue_T5_other               = FR_cue_T5 - FR_other;
    FR_cue_T5_other_Real          = FR_cue_T5_other(1,:);
    infoP_cue_T5                  = sum(repmat(FR_cue_T5_other_Real,numShuffles-1,1) >= FR_cue_T5_other(2:end,:),1)./(numShuffles-1);
    
    infoP_cue_T                   = []
    infoP_cue_T(1,:)              = infoP_cue;     
    infoP_cue_T(2,:)              = infoP_cue_T1;    
    infoP_cue_T(3,:)              = infoP_cue_T2;    
    infoP_cue_T(4,:)              = infoP_cue_T3;    
    infoP_cue_T(5,:)              = infoP_cue_T4;    
    infoP_cue_T(6,:)              = infoP_cue_T5;    
    
    Num_cue_095                   = sum(infoP_cue_T>=0.95,1)
    Num_cue_005                   = sum(infoP_cue_T<=0.05,1)
    
    tt=(infoP_cue_T>=0.95)
    
    infoP_cue_T                   = infoP_cue_T';
    
    %%-------------------------------------------------------------------------
    
    NumY_show                         = 1000


    for Cue_Trial                     = 1:size(infoP_cue_T,2)
    
    
        [infoP_cue_T_sort(:,Cue_Trial)  infoP_cue_T_Index(:,Cue_Trial)] = sortrows(infoP_cue_T(:,Cue_Trial),'descend');
        
        Cells_cueup(:,Cue_Trial)      = (infoP_cue_T_sort(:,Cue_Trial)>=0.95);
        Cells_cueup_Num(Cue_Trial)    = sum(double(Cells_cueup(:,Cue_Trial) ));
           
        Cells_cuedown(:,Cue_Trial)    = (infoP_cue_T_sort(:,Cue_Trial)<=0.05);
        Cells_cuedown_Num(Cue_Trial)  = sum(double(Cells_cuedown(:,Cue_Trial) ));
    
            
        Cells_cuebetween(:,Cue_Trial) = (infoP_cue_T_sort(:,Cue_Trial)<0.95)&(infoP_cue_T_sort(:,Cue_Trial)>0.05);
        
        
        C_trace_sort_good             = C_trace_unify_good(infoP_cue_T_Index(:,Cue_Trial),:);
        C_trace_sort_cueup            = C_trace_sort_good(Cells_cueup(:,Cue_Trial),:);
        C_trace_sort_cuedown          = C_trace_sort_good(Cells_cuedown(:,Cue_Trial),:);
        C_trace_sort_cuebetween       = C_trace_sort_good(Cells_cuebetween(:,Cue_Trial),:);
        
        
        fig                           = figure( 3000+Cue_Trial )
        fig.Position                  = [10 10 300 1.5*size(C_trace_sort_cueup,1)]; 

        %subplot(2,1,1);

        for i                         = 1:min(NumY_show, size(C_trace_sort_cueup,1))
            plot(C_trace_sort_cueup(i,:)*3+i);
            hold on;
        end

        axis([1, size(C_trace_sort_cueup,2), -2, min(NumY_show, size(C_trace_sort_cueup,1))+2]);
        

        for i  = 1:numel(miniscopecuebeginFrame)
            % 定义时间间隔
            t1 = miniscopecuebeginFrame(i); % 开始时间
            t2 = miniscopecueendFrame(i); % 结束时间
            % 添加阴影
            px = [t1, t2, t2, t1]; % 阴影的 x 坐标
            py = [-2, -2, min(NumY_show, size(C_trace_sort_cueup,1))+2, min(NumY_show, size(C_trace_sort_cueup,1))+2]; % 阴影的 y 坐标
            patch(px, py, 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none'); % 添加阴影
        
    %         line([t1 t1], [-2 min(NumY_show, size(C_trace_sort_cueup,1))+2], 'Color', 'r', 'LineStyle', '--'); % 添加开始时间的竖直线
    %         line([t2 t2], [-2 min(NumY_show, size(C_trace_sort_cueup,1))+2], 'Color', 'r', 'LineStyle', '--'); % 添加结束时间的竖直线
        end
        title('The first 1000 cue up cell trace')

       
        fig2          = figure(4000+Cue_Trial)
        fig2.Position = [500 10 300 1*size(C_trace_sort_cuedown,1)]; 

        %subplot(2,1,1);
        for i = 1:min(NumY_show, size(C_trace_sort_cuedown,1))
            plot(C_trace_sort_cuedown(i,:)*3+i);
            hold on;
        end

        axis([1, size(C_trace_sort_cuedown,2), -2, min(NumY_show, size(C_trace_sort_cuedown,1))+2]);
        
        for i  = 1:numel(miniscopecuebeginFrame)
            % 定义时间间隔
            t1 = miniscopecuebeginFrame(i); % 开始时间
            t2 = miniscopecueendFrame(i); % 结束时间
            % 添加阴影
            px = [t1, t2, t2, t1]; % 阴影的 x 坐标
            py = [-2, -2, min(NumY_show, size(C_trace_sort_cuedown,1))+2, min(NumY_show, size(C_trace_sort_cuedown,1))+2]; % 阴影的 y 坐标
            patch(px, py, 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none'); % 添加阴影
        
    %         line([t1 t1], [-2 min(NumY_show, size(C_trace_sort_cuedown,1))+2], 'Color', 'r', 'LineStyle', '--'); % 添加开始时间的竖直线
    %         line([t2 t2], [-2 min(NumY_show, size(C_trace_sort_cuedown,1))+2], 'Color', 'r', 'LineStyle', '--'); % 添加结束时间的竖直线
        end
    
      
    
        title('The first 1000 cue down cell trace')

         
    Cue_Trial
        
    end 


Color_cueup_trial        = [255 0     0;
                            255 165   0;
                            0   255   0;
                            0   0   255;
                            139 0   255]/255;    


for Cue_Trial                       = 2:size(infoP_cue_T,2)


    C_trace_sort_good             = C_trace_good(infoP_cue_T_Index(:,Cue_Trial),:);
    C_trace_sort_cueup            = C_trace_sort_good(Cells_cueup(:,Cue_Trial),:);
    C_trace_sort_cueup_Avg        = sum(C_trace_sort_cueup,1)/size(C_trace_sort_cueup,1);
    figure(4500)
    plot( C_trace_sort_cueup_Avg, 'Color',  Color_cueup_trial(Cue_Trial-1,:));
    hold on;

    h1 = fill(1:size(C_trace_sort_good,2), speed_COMET_smooth_7006_binary*20,'blue','EdgeColor','none');
    set(h1,'facealpha',.05)

%     b1=bar(speed_COMET_smooth_7006_binary*20,'FaceColor',[0 0 0])
%     b1.FaceAlpha = 0.3;
%     plot(-speed_COMET_smooth_7006  ,'Color', 'black')

    hold on;
    %         axis([1, size(C_trace_sort_cuedown,2), 0, 1]);
    
    for i=1:numel(miniscopecuebeginFrame)
        % 定义时间间隔
        t1 = miniscopecuebeginFrame(i); % 开始时间
        t2 = miniscopecueendFrame(i); % 结束时间
        % 添加阴影
        px = [t1, t2, t2, t1]; % 阴影的 x 坐标
        py = [0, 0,25, 25]; % 阴影的 y 坐标
        patch(px, py, 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none'); % 添加阴影

        %         line([t1 t1], [-2 min(NumY_show, size(C_trace_sort_cuedown,1))+2], 'Color', 'r', 'LineStyle', '--'); % 添加开始时间的竖直线
        %         line([t2 t2], [-2 min(NumY_show, size(C_trace_sort_cuedown,1))+2], 'Color', 'r', 'LineStyle', '--'); % 添加结束时间的竖直线
    end

    hold on
end


for Cue_Trial                 = 2:size(infoP_cue_T,2)


    C_trace_sort_good        = C_trace_good(infoP_cue_T_Index(:,Cue_Trial),:);
    C_trace_sort_cuedown     = C_trace_sort_good(Cells_cuedown(:,Cue_Trial),:);

    figure(4501)
    plot(  sum(C_trace_sort_cuedown,1)/size(C_trace_sort_cuedown,1), 'Color',  Color_cueup_trial(Cue_Trial-1,:));
    hold on;

    h1 = fill(1:size(C_trace_sort_good,2), speed_COMET_smooth_7006_binary*20,'blue','EdgeColor','none');
    set(h1,'facealpha',.05)

    %         plot(-speed_COMET_smooth_7006*20  ,'Color', 'black')
    hold on;

    %         axis([1, size(C_trace_sort_cuedown,2), 0, 1]);

    for i=1:numel(miniscopecuebeginFrame)
        % 定义时间间隔
        t1 = miniscopecuebeginFrame(i); % 开始时间
        t2 = miniscopecueendFrame(i); % 结束时间
        % 添加阴影
        px = [t1, t2, t2, t1]; % 阴影的 x 坐标
        py = [0, 0,25, 25]; % 阴影的 y 坐标
        patch(px, py, 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none'); % 添加阴影

        %         line([t1 t1], [-2 min(NumY_show, size(C_trace_sort_cuedown,1))+2], 'Color', 'r', 'LineStyle', '--'); % 添加开始时间的竖直线
        %         line([t2 t2], [-2 min(NumY_show, size(C_trace_sort_cuedown,1))+2], 'Color', 'r', 'LineStyle', '--'); % 添加结束时间的竖直线
    end

    hold on


end
%%===Analyze cue related neurons

%%=========================================================================


%% cueup spatial distrbution

Cell_cueup_id              = [];

Cell_cueup_id.T            = infoP_cue_T_Index(infoP_cue_T_sort(:,1)>=0.95,1);
Cell_cueup_id.T1           = infoP_cue_T_Index(infoP_cue_T_sort(:,2)>=0.95,2);
Cell_cueup_id.T2           = infoP_cue_T_Index(infoP_cue_T_sort(:,3)>=0.95,3);
Cell_cueup_id.T3           = infoP_cue_T_Index(infoP_cue_T_sort(:,4)>=0.95,4);
Cell_cueup_id.T4           = infoP_cue_T_Index(infoP_cue_T_sort(:,5)>=0.95,5);
Cell_cueup_id.T5           = infoP_cue_T_Index(infoP_cue_T_sort(:,6)>=0.95,6);

Cell_cueup_id.T_all        = {};
Cell_cueup_id.T_all{1}     = Cell_cueup_id.T; 
Cell_cueup_id.T_all{2}     = Cell_cueup_id.T1; 
Cell_cueup_id.T_all{3}     = Cell_cueup_id.T2; 
Cell_cueup_id.T_all{4}     = Cell_cueup_id.T3; 
Cell_cueup_id.T_all{5}     = Cell_cueup_id.T4; 
Cell_cueup_id.T_all{6}     = Cell_cueup_id.T5; 

Cell_cueup_id.T_index      = [];
P_thrd                     = 0.95;
Cell_cueup_id.T_index(:,1) = (infoP_cue_T(:,1) >=P_thrd);
Cell_cueup_id.T_index(:,2) = (infoP_cue_T(:,2) >=P_thrd);
Cell_cueup_id.T_index(:,3) = (infoP_cue_T(:,3) >=P_thrd);
Cell_cueup_id.T_index(:,4) = (infoP_cue_T(:,4) >=P_thrd);
Cell_cueup_id.T_index(:,5) = (infoP_cue_T(:,5) >=P_thrd);
Cell_cueup_id.T_index(:,6) = (infoP_cue_T(:,6) >=P_thrd);


for trialidx=1:6

    Cell_cueup_id.T_Num(trialidx)        = length(Cell_cueup_id.T_all{trialidx}) ;
    trialidx

end


Cell_cueup_id.T_index(:,7) = sum(Cell_cueup_id.T_index(:,5:6),2);
sum(  Cell_cueup_id.T_index(:,7)==2 )

find( Cell_cueup_id.T_index(:,7)==3 )

%%-------------Save the id of the cueup regulated cells--------------------

COMET_Correction.Cell_cueup_id  = Cell_cueup_id;
COMET_Correction.Cell_cueup_idx = ( infoP_cue_T>=0.95) ;

%%-------------Save the id of the cueup regulated cells--------------------

color                           = zeros(length(COMET_Correction.Cell_cueup_id.T) , 3);

unique_area_colormap            = []
for idx                         = 1:size(atlas_new.unique_area,2)

    unique_area_colormap(idx,:) = (atlas_new.unique_area{2,idx});
    unique_area_name{idx,:}     = (atlas_new.unique_area{1,idx});

end

%%---------------------Show the cueup cells from the first to last trial--

Cell_cueup_id.Numbrainregion   = zeros(size(Cell_cueup_id.T_all,2), length(atlas_new.cortical_neuron_id));

for trialidx = 1:size(Cell_cueup_id.T_all,2)
    
   color                             = []
   Cell_cueup_id_brainregion         = []

    for k    = 1:length(Cell_cueup_id.T_all{trialidx})
    
        brain_region_id              = COMET_Correction.CellRegion_id( Cell_cueup_id.T_all{trialidx}(k) );
        brain_region_name            = atlas_new.nameLR{brain_region_id};
        [~,unique_area_name_idx]     = ismember(unique_area_name,brain_region_name,'rows');
        color( k , : )               = unique_area_colormap( unique_area_name_idx>0 , :);
        Cell_cueup_id_brainregion(k) = brain_region_id;
        
    end
    
    COMET_Correction.Cell_cueup_id_brainregion{trialidx} = Cell_cueup_id_brainregion;
    
    A_listgood                       = (Cell_cueup_id.T_all{trialidx});

%     A_good_sparse            = A_neuron_sparse(:,A_neuron_good_idx+1);
%     A_good_sparse            = A_good_sparse(:,A_listgood);
%     A_good_2D_norm           = normalize(A_good_sparse,1,'range');
%     A_color                  = (A_good_2D_norm.*(A_good_2D_norm>0.3))*color;
%     A_color_2D               = reshape( A_color , 1944,1944,3 )*1;
%     
    color_mask               = [1,1,1];
%     
    pruned_atlas_mask_color  = reshape( atlas_new.pruned_atlas_mask_crop2COMET, 1944*1944 , 1 )*color_mask;
    pruned_atlas_mask_color  = reshape( pruned_atlas_mask_color , 1944,1944,3 )*1;

%     A_color_2D2              = (rot90(A_color_2D, -2))+pruned_atlas_mask_color;
    
%     figure(3000+trialidx)
%     imshow(A_color_2D2)%% not yet correct the coordinates of A in atlas!!!
%     
    N_Rcenter_r2             = COMET_Correction.N_Rcenter;
    
    figure( 5000+trialidx )
    imshow( pruned_atlas_mask_color )
    hold on;   
   
    for i=1:length(atlas_new.cortical_neuron_id)       
        
        cell_brainregion_idx                       = intersect( A_listgood , atlas_new.cortical_neuron_id{i} );
        scatter(N_Rcenter_r2(cell_brainregion_idx,1),N_Rcenter_r2(cell_brainregion_idx,2),8,atlas_new.unique_area{2,i} ,'filled')
        Cell_cueup_id.Numbrainregion( trialidx, i) = size(cell_brainregion_idx,1);
       
        axis equal

    end
    hold off;
 
%     % 用于存储图例条目的句柄
%     h = zeros(1, size(atlas_new.unique_area, 2));
%     % 遍历每种颜色，绘制带有对应颜色的圆圈
%     for i = 1:size(atlas_new.unique_area, 2)
%         colorRGB = atlas_new.unique_area{2, i}; % RGB 值
%         h(i) = plot(NaN,NaN,'o','MarkerEdgeColor',colorRGB,'MarkerFaceColor',colorRGB); % 使用 NaN 来创建图例条目
%     end
%     % 添加图例
%     lgd=legend(h, atlas_new.unique_area(1, :));
%     % 去除图例边框和背景
%     set(lgd, 'Box', 'off', 'Color', 'none');
%     lgd.TextColor = 'black';
%     % 隐藏坐标轴
%     axis off;
%     % exportgraphics(gcf, ['cueup','_P0.95.png'], 'Resolution', 500);

end


figure( 5100 )

color_mask               = [1,  1,   1];

Color_cueup_trial        = [255 0     0;
                            255 165   0;
                            0   255   0;
                            0   0   255;
                            139 0   255]/255;    

Marker_cueup_trial       = ['o'; '+'; '*' ;'x' ;'s'];

pruned_atlas_mask_color  = reshape( atlas_new.pruned_atlas_mask_crop2COMET, 1944*1944 , 1 )*color_mask;
pruned_atlas_mask_color  = reshape( pruned_atlas_mask_color , 1944,1944,3 )*1;

imshow( pruned_atlas_mask_color )
hold on;


for trialidx             = 2:size(Cell_cueup_id.T_all,2)
    
   color                     = []

   Cell_cueup_id_brainregion = []

    for k                            = 1:length(Cell_cueup_id.T_all{trialidx})
    
        brain_region_id              = COMET_Correction.CellRegion_id( Cell_cueup_id.T_all{trialidx}(k) );
        brain_region_name            = atlas_new.nameLR{brain_region_id};
        [~,unique_area_name_idx]     = ismember(unique_area_name,brain_region_name,'rows');
        color( k , : )               = unique_area_colormap( unique_area_name_idx>0 , :);
        Cell_cueup_id_brainregion(k) = brain_region_id;
        
    end
    
    COMET_Correction.Cell_cueup_id_brainregion{trialidx} = Cell_cueup_id_brainregion;
    
    A_listgood                       = (Cell_cueup_id.T_all{trialidx});
   
    for i=1:length(atlas_new.cortical_neuron_id)       
        
        cell_brainregion_idx                       = intersect( A_listgood , atlas_new.cortical_neuron_id{i} );
        scatter(N_Rcenter_r2(cell_brainregion_idx,1),N_Rcenter_r2(cell_brainregion_idx,2),16, Color_cueup_trial(trialidx-1,:),'Marker', Marker_cueup_trial(trialidx-1))
        Cell_cueup_id.Numbrainregion( trialidx, i) = size(cell_brainregion_idx,1);
       
        axis equal

    end
    hold on;
 

end





for i=1:length(atlas_new.cortical_neuron_id) 
    
    Cell_cueup_id.NumCellsall_brainregion(i) = size( atlas_new.cortical_neuron_id{1,i},2);
    
end

Cell_cueup_id.percent        = Cell_cueup_id.Numbrainregion ./ repmat( Cell_cueup_id.NumCellsall_brainregion, size(Cell_cueup_id.T_all,2), 1);




brainname                    = categorical(string(atlas_new.unique_area(1,:))');
brain_per                    = Cell_cueup_id.percent(2:6,:)' ;

figure(5500)
bar( brainname,   brain_per );

pruned_atlas_mask_color      = reshape( atlas_new.pruned_atlas_mask_crop2COMET, 1944*1944 , 1 )*color_mask;
pruned_atlas_mask_color      = reshape( pruned_atlas_mask_color , 1944,1944,3 )*1;
N_Rcenter_r2                 = COMET_Correction.N_Rcenter;

figure(6000)
imshow(pruned_atlas_mask_color)
hold on;

for trialidx = 1:size(Cell_cueup_id.T_all,2)

    A_listgood               = (Cell_cueup_id.T_all{trialidx});

    for i=1:length(atlas_new.cortical_neuron_id)

        cell_brainregion_idx = intersect( A_listgood , atlas_new.cortical_neuron_id{i} );
        scatter(N_Rcenter_r2(cell_brainregion_idx,1),N_Rcenter_r2(cell_brainregion_idx,2),trialidx*2,color(trialidx,:))
        axis equal    

    end
    hold on;
    trialidx

end

%--------------------------------------------------------------------------

COMET_Correction.Cell_cueup_id.T1


C_trace_sort_good        = C_trace_good(infoP_cue_T_Index(:,Cue_Trial),:);
        C_trace_sort_cuedown     = C_trace_sort_good(Cells_cuedown(:,Cue_Trial),:);



figure(6500)

plot(  sum( C_trace_good( Cell_cueup_id.T1,: ),1 )/size(Cell_cueup_id.T1,1),'Color', Color_cueup_trial(1,:)  );
hold on;
plot(  sum( C_trace_good( Cell_cueup_id.T2,: ),1 )/size(Cell_cueup_id.T2,1),'Color', Color_cueup_trial(2,:)   );
hold on;
plot(  sum( C_trace_good( Cell_cueup_id.T3,: ),1 )/size(Cell_cueup_id.T3,1),'Color', Color_cueup_trial(3,:)   );
hold on;
plot(  sum( C_trace_good( Cell_cueup_id.T4,: ),1 )/size(Cell_cueup_id.T4,1),'Color', Color_cueup_trial(4,:)   );
hold on;
plot(  sum( C_trace_good( Cell_cueup_id.T5,: ),1 )/size(Cell_cueup_id.T5,1),'Color', Color_cueup_trial(5,:)   );
hold on;

    h1 = fill(1:size(C_trace_sort_good,2), speed_COMET_smooth_7006_binary*20,'blue','EdgeColor','none');
    set(h1,'facealpha',.2)
hold on;
% plot(-speed_COMET_smooth_7006  ,'Color', 'black')
    hold on;
for i=1:numel(miniscopecuebeginFrame)
    % 定义时间间隔
    t1 = miniscopecuebeginFrame(i); % 开始时间
    t2 = miniscopecueendFrame(i); % 结束时间
    % 添加阴影
    px = [t1, t2, t2, t1]; % 阴影的 x 坐标
    py = [-2, -2, 12, 12]; % 阴影的 y 坐标
    patch(px, py, 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none'); % 添加阴影

    %         line([t1 t1], [-2 min(NumY_show, size(C_trace_sort_cuedown,1))+2], 'Color', 'r', 'LineStyle', '--'); % 添加开始时间的竖直线
    %         line([t2 t2], [-2 min(NumY_show, size(C_trace_sort_cuedown,1))+2], 'Color', 'r', 'LineStyle', '--'); % 添加结束时间的竖直线
end


figure(6501)

plot(  sum( C_trace_good( Cell_cueup_id.T,: ),1 )/size(Cell_cueup_id.T1,1),'black'  );
hold on;

    h1 = fill(1:size(C_trace_sort_good,2), speed_COMET_smooth_7006_binary*20,'blue','EdgeColor','none');
    set(h1,'facealpha',.2)
hold on;
% plot(-speed_COMET_smooth_7006  ,'Color', 'black')
hold on;
for i=1:numel(miniscopecuebeginFrame)
    % 定义时间间隔
    t1 = miniscopecuebeginFrame(i); % 开始时间
    t2 = miniscopecueendFrame(i); % 结束时间
    % 添加阴影
    px = [t1, t2, t2, t1]; % 阴影的 x 坐标
    py = [-2, -2, 12, 12]; % 阴影的 y 坐标
    patch(px, py, 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none'); % 添加阴影

    %         line([t1 t1], [-2 min(NumY_show, size(C_trace_sort_cuedown,1))+2], 'Color', 'r', 'LineStyle', '--'); % 添加开始时间的竖直线
    %         line([t2 t2], [-2 min(NumY_show, size(C_trace_sort_cuedown,1))+2], 'Color', 'r', 'LineStyle', '--'); % 添加结束时间的竖直线
end



%%---------------------Show the cueup cells from the first to last trial--
%%-------------------------------------------------------------------------
%%=========================================================================
%% cuedown spatial distrbution

Cell_cuedown_id              = [];
P_thrd                       = 0.05;

Cell_cuedown_id.T            = infoP_cue_T_Index( infoP_cue_T_sort(:,1)<=P_thrd,1 );
Cell_cuedown_id.T1           = infoP_cue_T_Index( infoP_cue_T_sort(:,2)<=P_thrd,2 );
Cell_cuedown_id.T2           = infoP_cue_T_Index( infoP_cue_T_sort(:,3)<=P_thrd,3 );
Cell_cuedown_id.T3           = infoP_cue_T_Index( infoP_cue_T_sort(:,4)<=P_thrd,4 );
Cell_cuedown_id.T4           = infoP_cue_T_Index( infoP_cue_T_sort(:,5)<=P_thrd,5 );
Cell_cuedown_id.T5           = infoP_cue_T_Index( infoP_cue_T_sort(:,6)<=P_thrd,6 );

Cell_cuedown_id.T_all        = {};
Cell_cuedown_id.T_all{1}     = Cell_cuedown_id.T; 
Cell_cuedown_id.T_all{2}     = Cell_cuedown_id.T1; 
Cell_cuedown_id.T_all{3}     = Cell_cuedown_id.T2; 
Cell_cuedown_id.T_all{4}     = Cell_cuedown_id.T3; 
Cell_cuedown_id.T_all{5}     = Cell_cuedown_id.T4; 
Cell_cuedown_id.T_all{6}     = Cell_cuedown_id.T5; 
 
Cell_cuedown_id.T_index      = [];

Cell_cuedown_id.T_index(:,1) = (infoP_cue_T(:,1) <=P_thrd);
Cell_cuedown_id.T_index(:,2) = (infoP_cue_T(:,2) <=P_thrd);
Cell_cuedown_id.T_index(:,3) = (infoP_cue_T(:,3) <=P_thrd);
Cell_cuedown_id.T_index(:,4) = (infoP_cue_T(:,4) <=P_thrd);
Cell_cuedown_id.T_index(:,5) = (infoP_cue_T(:,5) <=P_thrd);
Cell_cuedown_id.T_index(:,6) = (infoP_cue_T(:,6) <=P_thrd);


for trialidx=1:6

    Cell_cuedown_id.T_Num(trialidx)        = length(Cell_cuedown_id.T_all{trialidx}) ;
    trialidx

end


Cell_cuedown_id.T_index(:,7) = sum(Cell_cuedown_id.T_index(:,5:6),2);
sum(  Cell_cuedown_id.T_index(:,7)==2 )

%%-------------Save the id of the cueup regulated cells--------------------

COMET_Correction.Cell_cuedown_id  = Cell_cuedown_id;
COMET_Correction.Cell_cuedown_idx = ( infoP_cue_T<=P_thrd) ;

%%-------------Save the id of the cueup regulated cells--------------------
color                             = zeros(length(COMET_Correction.Cell_cuedown_id.T) , 3);
%%---------------------Show the cueup cells from the first to last trial--

for trialidx = 1:size(Cell_cuedown_id.T_all,2)
    
    color                              = []
    Cell_cuedown_id_brainregion        = []

    for k=1:length(Cell_cuedown_id.T_all{trialidx})

        brain_region_id                = COMET_Correction.CellRegion_id( Cell_cuedown_id.T_all{trialidx}(k) );
        brain_region_name              = atlas_new.nameLR{brain_region_id};
        [~,unique_area_name_idx]       = ismember(unique_area_name,brain_region_name,'rows');
        color( k , : )                 = unique_area_colormap( unique_area_name_idx>0 , :);
        Cell_cuedown_id_brainregion(k) =  brain_region_id;

    end

    COMET_Correction.Cell_cuedown_id_brainregion{trialidx}=Cell_cuedown_id_brainregion;

    A_listgood                         = (Cell_cuedown_id.T_all{trialidx});

    %     A_good_sparse            = A_neuron_sparse(:,A_neuron_good_idx+1);
    %     A_good_sparse            = A_good_sparse(:,A_listgood);
    %     A_good_2D_norm           = normalize(A_good_sparse,1,'range');
    %     A_color                  = (A_good_2D_norm.*(A_good_2D_norm>0.3))*color;
    %     A_color_2D               = reshape( A_color , 1944,1944,3 )*1;

    color_mask                         = [1,1,1];

    %     pruned_atlas_mask_color  = reshape( atlas_new.pruned_atlas_mask_crop2COMET, 1944*1944 , 1 )*color_mask;
    %     pruned_atlas_mask_color  = reshape( pruned_atlas_mask_color , 1944,1944,3 )*1;
    %     A_color_2D2              = (rot90(A_color_2D, -2)) + pruned_atlas_mask_color;
    %
    %     figure(3000+trialidx)
    %     imshow(A_color_2D2)%% not yet correct the coordinates of A in atlas!!!
    %

    N_Rcenter_r2                   = COMET_Correction.N_Rcenter;

    figure( 7000 + trialidx )
    imshow(  pruned_atlas_mask_color  )
    hold on;
    for i=1:length(atlas_new.cortical_neuron_id)

        cell_brainregion_idx = intersect( A_listgood , atlas_new.cortical_neuron_id{i} );
        scatter(N_Rcenter_r2(cell_brainregion_idx,1),N_Rcenter_r2(cell_brainregion_idx,2),8,atlas_new.unique_area{2,i} ,'filled')

        Cell_cuedown_id.Numbrainregion( trialidx, i) = size(cell_brainregion_idx,1);

        axis equal

    end
    hold off;


    %     % 用于存储图例条目的句柄
    %     h = zeros(1, size(atlas_new.unique_area, 2));
    %     % 遍历每种颜色，绘制带有对应颜色的圆圈
    %     for i = 1:size(atlas_new.unique_area, 2)
    %         colorRGB = atlas_new.unique_area{2, i}; % RGB 值
    %         h(i) = plot(NaN,NaN,'o','MarkerEdgeColor',colorRGB,'MarkerFaceColor',colorRGB); % 使用 NaN 来创建图例条目
    %     end
    %     % 添加图例
    %     lgd=legend(h, atlas_new.unique_area(1, :));
    %     % 去除图例边框和背景
    %     set(lgd, 'Box', 'off', 'Color', 'none');
    %     lgd.TextColor = 'black';
    %     % 隐藏坐标轴
    %     axis off;
    %     % exportgraphics(gcf, ['cueup','_P0.95.png'], 'Resolution', 500);

end



Color_cuedown_trial        = [255 0     0;
                              255 165   0;
                              0   255   0;
                              0   0   255;
                              139 0   255]/255;    

Marker_cuedown_trial       = ['o'; '+'; '*' ;'x' ;'s'];

    figure( 7100)
    imshow(  pruned_atlas_mask_color  )
hold on
for trialidx = 2:size(Cell_cuedown_id.T_all,2)
    

    A_listgood                     = (Cell_cuedown_id.T_all{trialidx});
  
    N_Rcenter_r2                   = COMET_Correction.N_Rcenter;


    for i=1:length(atlas_new.cortical_neuron_id)

        cell_brainregion_idx = intersect( A_listgood , atlas_new.cortical_neuron_id{i} );
        scatter(N_Rcenter_r2(cell_brainregion_idx,1),N_Rcenter_r2(cell_brainregion_idx,2),16, Color_cuedown_trial(trialidx-1,:),'Marker', Marker_cuedown_trial(trialidx-1))

          axis equal

    end
    hold on;

   
end



for i=1:length(atlas_new.cortical_neuron_id) 
    
    Cell_cuedown_id.NumCellsall_brainregion(i) = size( atlas_new.cortical_neuron_id{1,i},2);
    
end

Cell_cuedown_id.percent        = Cell_cuedown_id.Numbrainregion ./ repmat( Cell_cuedown_id.NumCellsall_brainregion, size(Cell_cuedown_id.T_all,2), 1);
brainname                      = categorical(string(atlas_new.unique_area(1,:))');
brain_cuedown_per              = Cell_cuedown_id.percent(2:6,:)' ;

figure(7500)
bar( brainname,   brain_cuedown_per );

%--------------------------------------------------------------------------

figure(8000)

plot(  sum( C_trace_good( Cell_cuedown_id.T1,: ),1 )/size(Cell_cuedown_id.T1,1),'Color',Color_cuedown_trial(1,:)  );
hold on;
plot(  sum( C_trace_good( Cell_cuedown_id.T2,: ),1 )/size(Cell_cuedown_id.T2,1),'Color',Color_cuedown_trial(2,:)  );
hold on;
plot(  sum( C_trace_good( Cell_cuedown_id.T3,: ),1 )/size(Cell_cuedown_id.T3,1),'Color',Color_cuedown_trial(3,:)  );
hold on;
plot(  sum( C_trace_good( Cell_cuedown_id.T4,: ),1 )/size(Cell_cuedown_id.T4,1),'Color',Color_cuedown_trial(4,:)  );
hold on;
plot(  sum( C_trace_good( Cell_cuedown_id.T5,: ),1 )/size(Cell_cuedown_id.T5,1),'Color',Color_cuedown_trial(5,:)  );
hold on;

    h1 = fill(1:size(C_trace_sort_good,2), speed_COMET_smooth_7006_binary*20,'blue','EdgeColor','none');
    set(h1,'facealpha',.3)
hold on;
% plot(-speed_COMET_smooth_7006  ,'Color', 'black')
hold on;
for i=1:numel(miniscopecuebeginFrame)
    % 定义时间间隔
    t1 = miniscopecuebeginFrame(i); % 开始时间
    t2 = miniscopecueendFrame(i); % 结束时间
    % 添加阴影
    px = [t1, t2, t2, t1]; % 阴影的 x 坐标
    py = [-2, -2, 12, 12]; % 阴影的 y 坐标
    patch(px, py, 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none'); % 添加阴影

    %         line([t1 t1], [-2 min(NumY_show, size(C_trace_sort_cuedown,1))+2], 'Color', 'r', 'LineStyle', '--'); % 添加开始时间的竖直线
    %         line([t2 t2], [-2 min(NumY_show, size(C_trace_sort_cuedown,1))+2], 'Color', 'r', 'LineStyle', '--'); % 添加结束时间的竖直线
end




figure(8001)

plot(  sum( C_trace_good( Cell_cuedown_id.T,: ),1 )/size(Cell_cuedown_id.T,1),'black'  );
hold on;

h1 = fill(1:size(C_trace_sort_good,2), speed_COMET_smooth_7006_binary*20,'blue','EdgeColor','none');
set(h1,'facealpha',.3)
hold on;
%   plot(-speed_COMET_smooth  ,'Color', 'black')
hold on;
for i=1:numel(miniscopecuebeginFrame)
    % 定义时间间隔
    t1 = miniscopecuebeginFrame(i); % 开始时间
    t2 = miniscopecueendFrame(i); % 结束时间
    % 添加阴影
    px = [t1, t2, t2, t1]; % 阴影的 x 坐标
    py = [-2, -2, 12, 12]; % 阴影的 y 坐标
    patch(px, py, 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none'); % 添加阴影

    %         line([t1 t1], [-2 min(NumY_show, size(C_trace_sort_cuedown,1))+2], 'Color', 'r', 'LineStyle', '--'); % 添加开始时间的竖直线
    %         line([t2 t2], [-2 min(NumY_show, size(C_trace_sort_cuedown,1))+2], 'Color', 'r', 'LineStyle', '--'); % 添加结束时间的竖直线
end





Color_cuedown_trial                    = [255 0     0;
                                          255 165   0;
                                          0   255   0;
                                          0   0   255;
                                          139 0   255]/255;    

Marker_cuedown_trial                   = ['o'; '+'; '*' ;'x' ;'s'];


    figure( 9100)
    imshow(  pruned_atlas_mask_color  )
hold on
for trialidx = 2:size(Cell_cuedown_id.T_all,2)
    
    color                              = []
    Cell_cuedown_id_brainregion        = []

    for k=1:length(Cell_cuedown_id.T_all{trialidx})

        brain_region_id                = COMET_Correction.CellRegion_id( Cell_cuedown_id.T_all{trialidx}(k) );
        brain_region_name              = atlas_new.nameLR{brain_region_id};
        [~,unique_area_name_idx]       = ismember(unique_area_name,brain_region_name,'rows');
        color( k , : )                 = unique_area_colormap( unique_area_name_idx>0 , :);
        Cell_cuedown_id_brainregion(k) =  brain_region_id;

    end

    COMET_Correction.Cell_cuedown_id_brainregion{trialidx}=Cell_cuedown_id_brainregion;

    A_listgood                         = (Cell_cuedown_id.T_all{trialidx});

  
    color_mask                         = [1,1,1];

  
    N_Rcenter_r2                       = COMET_Correction.N_Rcenter;


    for i=1:length(atlas_new.cortical_neuron_id)

        cell_brainregion_idx = intersect( A_listgood , atlas_new.cortical_neuron_id{i} );
        scatter(N_Rcenter_r2(cell_brainregion_idx,1),N_Rcenter_r2(cell_brainregion_idx,2),16, Color_cuedown_trial(trialidx-1,:),'Marker', Marker_cuedown_trial(trialidx-1))

        Cell_cuedown_id.Numbrainregion( trialidx, i) = size(cell_brainregion_idx,1);

        axis equal

    end

    hold on;

   
end

%%=========================================================================
%%=========================================================================

%%===Analyze shock related neurons

Trace4Shock                       = C_trace_unify_good(:,COMET_Correction.ShockplusOther>0);
Frameidx_ShockplusOther           = COMET_Correction.ShockplusOther(COMET_Correction.ShockplusOther~=0);
figure(700)
plot(Frameidx_ShockplusOther)

ns                                = size(Trace4Shock,1);
numShuffles                       = 500;
numShuffles                       = numShuffles+1;
FR_Shock                          = nan(numShuffles,ns);
FR_Shock_T1                       = nan(numShuffles,ns);
FR_Shock_T2                       = nan(numShuffles,ns);
FR_Shock_T3                       = nan(numShuffles,ns);
FR_Shock_T4                       = nan(numShuffles,ns);
FR_Shock_T5                       = nan(numShuffles,ns);
FR_other                          = nan(numShuffles,ns);
NoFrames                          = size(Trace4Shock,2);
shift                             = floor( rand(numShuffles,1) *NoFrames);
shift(1,1)                        = 0;
indexp                            = 1:NoFrames;

T_spk_diff                        = [0;diff(timestamps_)];

T_spk_diff_Shockother             = T_spk_diff(COMET_Correction.ShockplusOther>0);
T_spk_diff_Shock                  = T_spk_diff(COMET_Correction.shockframe>0  );
T_spk_diff_Other                  = T_spk_diff(COMET_Correction.Other       );

T_spk_diff_shock_sum              = sum(T_spk_diff_Shock);
T_spk_diff_Other_sum              = sum(T_spk_diff_Other);

parfor idxnumShf                  = 1:numShuffles

    indexs2                       = mod(indexp+shift(idxnumShf), NoFrames);
    indexs2(indexs2==0)           = NoFrames;
    Trace4Shock_shf               = Trace4Shock(:,indexs2);

    Trace4Shock_shf_shock         = Trace4Shock_shf(:,Frameidx_ShockplusOther>=2);
    Trace4Shock_shf_shock_T1      = Trace4Shock_shf(:,Frameidx_ShockplusOther==2);
    Trace4Shock_shf_shock_T2      = Trace4Shock_shf(:,Frameidx_ShockplusOther==4);
    Trace4Shock_shf_shock_T3      = Trace4Shock_shf(:,Frameidx_ShockplusOther==6);
    Trace4Shock_shf_shock_T4      = Trace4Shock_shf(:,Frameidx_ShockplusOther==8);
    Trace4Shock_shf_shock_T5      = Trace4Shock_shf(:,Frameidx_ShockplusOther==10);
    Trace4Shock_shf_other         = Trace4Shock_shf(:,Frameidx_ShockplusOther==1);
    Trace4Shock_shf_shock_sum     = sum(Trace4Shock_shf_shock,2);
    Trace4Shock_shf_shock_T1_sum  = sum(Trace4Shock_shf_shock_T1,2);
    Trace4Shock_shf_shock_T2_sum  = sum(Trace4Shock_shf_shock_T2,2);
    Trace4Shock_shf_shock_T3_sum  = sum(Trace4Shock_shf_shock_T3,2);
    Trace4Shock_shf_shock_T4_sum  = sum(Trace4Shock_shf_shock_T4,2);
    Trace4Shock_shf_shock_T5_sum  = sum(Trace4Shock_shf_shock_T5,2);
    Trace4Shock_shf_other_sum     = sum(Trace4Shock_shf_other,2);

    FR_Shock(idxnumShf,:)         = Trace4Shock_shf_shock_sum/T_spk_diff_shock_sum;
    FR_Shock_T1(idxnumShf,:)      = Trace4Shock_shf_shock_T1_sum/T_spk_diff_shock_sum * 5;
    FR_Shock_T2(idxnumShf,:)      = Trace4Shock_shf_shock_T2_sum/T_spk_diff_shock_sum * 5;
    FR_Shock_T3(idxnumShf,:)      = Trace4Shock_shf_shock_T3_sum/T_spk_diff_shock_sum * 5;
    FR_Shock_T4(idxnumShf,:)      = Trace4Shock_shf_shock_T4_sum/T_spk_diff_shock_sum * 5;
    FR_Shock_T5(idxnumShf,:)      = Trace4Shock_shf_shock_T5_sum/T_spk_diff_shock_sum * 5;
    FR_othernoshock(idxnumShf,:)  = Trace4Shock_shf_other_sum/T_spk_diff_Other_sum;

    idxnumShf
end

FR_Shock_other                    = FR_Shock - FR_othernoshock;
FR_Shock_other_Real               = FR_Shock_other(1,:);
infoP_Shock                       = sum(repmat(FR_Shock_other_Real,numShuffles-1,1) >= FR_Shock_other(2:end,:),1)./(numShuffles-1);

FR_Shock_T1_other                 = FR_Shock_T1 - FR_othernoshock;
FR_Shock_T1_other_Real            = FR_Shock_T1_other(1,:);
infoP_Shock_T1                    = sum(repmat(FR_Shock_T1_other_Real,numShuffles-1,1) >= FR_Shock_T1_other(2:end,:),1)./(numShuffles-1);

FR_Shock_T2_other                 = FR_Shock_T2 - FR_othernoshock;
FR_Shock_T2_other_Real            = FR_Shock_T2_other(1,:);
infoP_Shock_T2                    = sum(repmat(FR_Shock_T2_other_Real,numShuffles-1,1) >= FR_Shock_T2_other(2:end,:),1)./(numShuffles-1);

FR_Shock_T3_other                 = FR_Shock_T3 - FR_othernoshock;
FR_Shock_T3_other_Real            = FR_Shock_T3_other(1,:);
infoP_Shock_T3                    = sum(repmat(FR_Shock_T3_other_Real,numShuffles-1,1) >= FR_Shock_T3_other(2:end,:),1)./(numShuffles-1);

FR_Shock_T4_other                 = FR_Shock_T4 - FR_othernoshock;
FR_Shock_T4_other_Real            = FR_Shock_T4_other(1,:);
infoP_Shock_T4                    = sum(repmat(FR_Shock_T4_other_Real,numShuffles-1,1) >= FR_Shock_T4_other(2:end,:),1)./(numShuffles-1);

FR_Shock_T5_other                 = FR_Shock_T5 - FR_othernoshock;
FR_Shock_T5_other_Real            = FR_Shock_T5_other(1,:);
infoP_Shock_T5                    = sum(repmat(FR_Shock_T5_other_Real,numShuffles-1,1) >= FR_Shock_T5_other(2:end,:),1)./(numShuffles-1);



infoP_Shock_T                     = []
infoP_Shock_T(1,:)                = infoP_Shock;     
infoP_Shock_T(2,:)                = infoP_Shock_T1;    
infoP_Shock_T(3,:)                = infoP_Shock_T2;    
infoP_Shock_T(4,:)                = infoP_Shock_T3;    
infoP_Shock_T(5,:)                = infoP_Shock_T4;    
infoP_Shock_T(6,:)                = infoP_Shock_T5;    

Num_Shock_095                     = sum( infoP_Shock_T>=0.95,1 )
Num_Shock_005                     = sum( infoP_Shock_T<=0.05,1 )

infoP_Shock_T                     = infoP_Shock_T';

%%-------------------------------------------------------------------------

NumY_show                         = 1000

infoP_Shock_T_sort                    = []
infoP_Shock_T_Index                   = []
Cells_Shockup                         = []



for Shock_Trial                       = 1:size(infoP_Shock_T,2)


    [infoP_Shock_T_sort(:,Shock_Trial)  infoP_Shock_T_Index(:,Shock_Trial)]  = sortrows(infoP_Shock_T(:, Shock_Trial),'descend');
    
    Cells_Shockup(:,Shock_Trial)      = (infoP_Shock_T_sort(:,Shock_Trial)>=0.95);
    Cells_Shockup_Num(Shock_Trial)    = sum(double(Cells_Shockup(:,Shock_Trial) ));
       
    Cells_Shockdown(:,Shock_Trial)    = (infoP_Shock_T_sort(:,Shock_Trial)<=0.05);
    Cells_Shockdown_Num(Shock_Trial)  = sum(double(Cells_Shockdown(:,Shock_Trial) ));

        
    Cells_Shockbetween(:,Shock_Trial) = (infoP_Shock_T_sort(:,Shock_Trial)<0.95)&(infoP_Shock_T_sort(:,Shock_Trial)>0.05);
    
    
    C_trace_sort_good                 = C_trace_unify_good(infoP_Shock_T_Index(:,Shock_Trial),:);
    C_trace_sort_Shockup              = C_trace_sort_good(logical(Cells_Shockup(:,Shock_Trial)),:);
    C_trace_sort_Shockdown            = C_trace_sort_good(Cells_Shockdown(:,Shock_Trial),:);
    C_trace_sort_Shockbetween         = C_trace_sort_good(Cells_Shockbetween(:,Shock_Trial),:);
    
    
    fig = figure( 8000+Shock_Trial )
    fig.Position = [10 10 300 1.5*size(C_trace_sort_Shockup,1)]; 
    %subplot(2,1,1);
    for i = 1:min(NumY_show, size(C_trace_sort_Shockup,1))
        plot(C_trace_sort_Shockup(i,:)*3+i,'k');
        hold on;
    end

    axis([1, size(C_trace_sort_Shockup,2), -2, min(NumY_show, size(C_trace_sort_Shockup,1))+2]);
    
    for i=1:numel(miniscopeshockbeginFrame)
        % 定义时间间隔
        t1 = miniscopeshockbeginFrame(i); % 开始时间
        t2 = miniscopeshockendFrame(i); % 结束时间
        % 添加阴影
        px = [t1, t2, t2, t1]; % 阴影的 x 坐标
        py = [-2, -2, min(NumY_show, size(C_trace_sort_Shockup,1))+2, min(NumY_show, size(C_trace_sort_Shockup,1))+2]; % 阴影的 y 坐标
        patch(px, py, 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none'); % 添加阴影
    
%         line([t1 t1], [-2 min(NumY_show, size(C_trace_sort_Shockup,1))+2], 'Color', 'r', 'LineStyle', '--'); % 添加开始时间的竖直线
%         line([t2 t2], [-2 min(NumY_show, size(C_trace_sort_Shockup,1))+2], 'Color', 'r', 'LineStyle', '--'); % 添加结束时间的竖直线
    end
    title('The first 1000 Shock up cell trace')


    fig2 = figure(9000+Shock_Trial)
    fig2.Position = [500 10 300 1*size(C_trace_sort_Shockdown,1)]; 
    %subplot(2,1,1);
    for i = 1:min(NumY_show, size(C_trace_sort_Shockdown,1))
        plot(C_trace_sort_Shockdown(i,:)*3+i,'k');
        hold on;
    end
    axis([1, size(C_trace_sort_Shockdown,2), -2, min(NumY_show, size(C_trace_sort_Shockdown,1))+2]);
    
    for i=1:numel(miniscopeshockbeginFrame)
        % 定义时间间隔
        t1 = miniscopeshockbeginFrame(i); % 开始时间
        t2 = miniscopeshockendFrame(i); % 结束时间
        % 添加阴影
        px = [t1, t2, t2, t1]; % 阴影的 x 坐标
        py = [-2, -2, min(NumY_show, size(C_trace_sort_Shockdown,1))+2, min(NumY_show, size(C_trace_sort_Shockdown,1))+2]; % 阴影的 y 坐标
        patch(px, py, 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none'); % 添加阴影
    
%         line([t1 t1], [-2 min(NumY_show, size(C_trace_sort_Shockdown,1))+2], 'Color', 'r', 'LineStyle', '--'); % 添加开始时间的竖直线
%         line([t2 t2], [-2 min(NumY_show, size(C_trace_sort_Shockdown,1))+2], 'Color', 'r', 'LineStyle', '--'); % 添加结束时间的竖直线
    end

    title('The first 1000 Shock down cell trace')

Shock_Trial
    
end 


%%=========================================================================
Cell_Shockup_id                   = [];
P_thrd                            = 0.95;

 

Cell_Shockup_id.T                 = infoP_Shock_T_Index( infoP_Shock_T_sort(:,1)>=P_thrd,1 );
Cell_Shockup_id.T1                = infoP_Shock_T_Index( infoP_Shock_T_sort(:,2)>=P_thrd,2 );
Cell_Shockup_id.T2                = infoP_Shock_T_Index( infoP_Shock_T_sort(:,3)>=P_thrd,3 );
Cell_Shockup_id.T3                = infoP_Shock_T_Index( infoP_Shock_T_sort(:,4)>=P_thrd,4 );
Cell_Shockup_id.T4                = infoP_Shock_T_Index( infoP_Shock_T_sort(:,5)>=P_thrd,5 );
Cell_Shockup_id.T5                = infoP_Shock_T_Index( infoP_Shock_T_sort(:,6)>=P_thrd,6 );

Cell_Shockup_id.T_all             = {};
Cell_Shockup_id.T_all{1}          = Cell_Shockup_id.T; 
Cell_Shockup_id.T_all{2}          = Cell_Shockup_id.T1; 
Cell_Shockup_id.T_all{3}          = Cell_Shockup_id.T2; 
Cell_Shockup_id.T_all{4}          = Cell_Shockup_id.T3; 
Cell_Shockup_id.T_all{5}          = Cell_Shockup_id.T4; 
Cell_Shockup_id.T_all{6}          = Cell_Shockup_id.T5; 

Cell_Shockup_id.T_index           = [];

Cell_Shockup_id.T_index(:,1)      = (infoP_Shock_T(:,1) >=P_thrd);
Cell_Shockup_id.T_index(:,2)      = (infoP_Shock_T(:,2) >=P_thrd);
Cell_Shockup_id.T_index(:,3)      = (infoP_Shock_T(:,3) >=P_thrd);
Cell_Shockup_id.T_index(:,4)      = (infoP_Shock_T(:,4) >=P_thrd);
Cell_Shockup_id.T_index(:,5)      = (infoP_Shock_T(:,5) >=P_thrd);
Cell_Shockup_id.T_index(:,6)      = (infoP_Shock_T(:,6) >=P_thrd);

for trialidx=1:6

    Cell_Shockup_id.T_Num(trialidx)        = length(Cell_Shockup_id.T_all{trialidx}) ;
    trialidx

end

for trialidx                      = 2:5

    Cell_Shockup_id.T_index(:,7)  = sum(Cell_Shockup_id.T_index(:,trialidx:trialidx+1),2);
    Cell_Shockup_id.T_index_Num(trialidx) = sum(  Cell_Shockup_id.T_index(:,7)==2 );
    trialidx

end

find( Cell_Shockup_id.T_index(:,7)==3 )


%%-------------Save the id of the Shockup regulated cells--------------------

COMET_Correction.Cell_Shockup_id  = Cell_Shockup_id;
COMET_Correction.Cell_Shockup_idx = ( infoP_Shock_T>=0.95) ;

%%-------------Save the id of the Shockup regulated cells--------------------
color                             = zeros(length(COMET_Correction.Cell_Shockup_id.T) , 3);
%%---------------------Show the Shockup cells from the first to last trial--

for trialidx = 1:size(Cell_Shockup_id.T_all,2)
    
   color                               = []
   Cell_Shockup_id_brainregion         = []

    for k=1:length(Cell_Shockup_id.T_all{trialidx})
    
        brain_region_id                = COMET_Correction.CellRegion_id( Cell_Shockup_id.T_all{trialidx}(k) );
        brain_region_name              = atlas_new.nameLR{brain_region_id};
        [~,unique_area_name_idx]       = ismember(unique_area_name,brain_region_name,'rows');
        color( k , : )                 = unique_area_colormap( unique_area_name_idx>0 , :);
        Cell_Shockup_id_brainregion(k) = brain_region_id;
        
    end
    
    COMET_Correction.Cell_Shockup_id_brainregion{trialidx}=Cell_Shockup_id_brainregion;
    
    A_listgood               = ( Cell_Shockup_id.T_all{trialidx} );
%     A_good_sparse            = A_neuron_sparse(:,A_neuron_good_idx+1);
%     A_good_sparse            = A_good_sparse(:,A_listgood);
%     A_good_2D_norm           = normalize(A_good_sparse,1,'range');
%     A_color                  = (A_good_2D_norm.*(A_good_2D_norm>0.3))*color;
%     A_color_2D               = reshape( A_color , 1944,1944,3 )*1;
    
    color_mask               = [1,1,1];
    
%     pruned_atlas_mask_color  = reshape( atlas_new.pruned_atlas_mask_crop2COMET, 1944*1944 , 1 )*color_mask;
%     pruned_atlas_mask_color  = reshape( pruned_atlas_mask_color , 1944,1944,3 )*1;
%     A_color_2D2              = (rot90(A_color_2D, -2))+pruned_atlas_mask_color;
%     
%     figure(3000+trialidx)
%     imshow(A_color_2D2)%% not yet correct the coordinates of A in atlas!!!
%     
    N_Rcenter_r2             = COMET_Correction.N_Rcenter;
    
    figure(9500+trialidx)
    imshow(pruned_atlas_mask_color)
    hold on;
    for i=1:length(atlas_new.cortical_neuron_id)       
        
        cell_brainregion_idx = intersect( A_listgood , atlas_new.cortical_neuron_id{i} );
        scatter(N_Rcenter_r2(cell_brainregion_idx,1),N_Rcenter_r2(cell_brainregion_idx,2),8,atlas_new.unique_area{2,i} ,'filled')

        Cell_Shockup_id.Numbrainregion( trialidx, i) = size(cell_brainregion_idx,1);

        axis equal

    end
    hold off;

 

%     % 用于存储图例条目的句柄
%     h = zeros(1, size(atlas_new.unique_area, 2));
%     % 遍历每种颜色，绘制带有对应颜色的圆圈
%     for i = 1:size(atlas_new.unique_area, 2)
%         colorRGB = atlas_new.unique_area{2, i}; % RGB 值
%         h(i) = plot(NaN,NaN,'o','MarkerEdgeColor',colorRGB,'MarkerFaceColor',colorRGB); % 使用 NaN 来创建图例条目
%     end
%     % 添加图例
%     lgd=legend(h, atlas_new.unique_area(1, :));
%     % 去除图例边框和背景
%     set(lgd, 'Box', 'off', 'Color', 'none');
%     lgd.TextColor = 'black';
%     % 隐藏坐标轴
%     axis off;
%     % exportgraphics(gcf, ['Shockup','_P0.95.png'], 'Resolution', 500);

end





for i=1:length(atlas_new.cortical_neuron_id) 
    
    Cell_Shockup_id.NumCellsall_brainregion(i) = size( atlas_new.cortical_neuron_id{1,i},2);
    
end

Cell_Shockup_id.percent        = Cell_Shockup_id.Numbrainregion ./ repmat( Cell_Shockup_id.NumCellsall_brainregion, size(Cell_Shockup_id.T_all,2), 1);
%--------------------------------------------------------------------------
Color_shockup_trial        = [255 0     0;
                              255 165   0;
                              0   255   0;
                              0   0   255;
                              139 0   255]/255;    

Marker_shockup_trial       = ['o'; '+'; '*' ;'x' ;'s'];




figure(8600)

plot(  sum( C_trace_good( Cell_Shockup_id.T1,: ),1 )/size(Cell_Shockup_id.T1,1),'Color', Color_shockup_trial(1,:)  );
hold on;
plot(  sum( C_trace_good( Cell_Shockup_id.T2,: ),1 )/size(Cell_Shockup_id.T2,1),'Color', Color_shockup_trial(2,:)  );
hold on;
plot(  sum( C_trace_good( Cell_Shockup_id.T3,: ),1 )/size(Cell_Shockup_id.T3,1),'Color', Color_shockup_trial(3,:)  );
hold on;
plot(  sum( C_trace_good( Cell_Shockup_id.T4,: ),1 )/size(Cell_Shockup_id.T4,1),'Color', Color_shockup_trial(4,:)  );
hold on;
plot(  sum( C_trace_good( Cell_Shockup_id.T5,: ),1 )/size(Cell_Shockup_id.T5,1),'Color', Color_shockup_trial(5,:)  );

hold on;
h1 = fill(1:size(C_trace_sort_good,2), speed_COMET_smooth_7006_binary*20,'blue','EdgeColor','none'             );
set(h1,'facealpha',.3)
hold on;
%     plot(-speed_COMET_smooth  ,'Color', 'black')
hold on;

for i=1:numel(miniscopecuebeginFrame)
    % 定义时间间隔
    t1 = miniscopeshockbeginFrame(i); % 开始时间
    t2 = miniscopeshockendFrame(i); % 结束时间
    % 添加阴影
    px = [t1, t2, t2, t1]; % 阴影的 x 坐标
    py = [-2, -2, 12, 12]; % 阴影的 y 坐标
    patch(px, py, 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none'); % 添加阴影

    %         line([t1 t1], [-2 min(NumY_show, size(C_trace_sort_cuedown,1))+2], 'Color', 'r', 'LineStyle', '--'); % 添加开始时间的竖直线
    %         line([t2 t2], [-2 min(NumY_show, size(C_trace_sort_cuedown,1))+2], 'Color', 'r', 'LineStyle', '--'); % 添加结束时间的竖直线
end


Cell_Shockdown_id                 = [];
P_thrd                            = 0.05;

Cell_Shockdown_id.T               = infoP_Shock_T_Index( infoP_Shock_T_sort(:,1)<=P_thrd,1 );
Cell_Shockdown_id.T1              = infoP_Shock_T_Index( infoP_Shock_T_sort(:,2)<=P_thrd,2 );
Cell_Shockdown_id.T2              = infoP_Shock_T_Index( infoP_Shock_T_sort(:,3)<=P_thrd,3 );
Cell_Shockdown_id.T3              = infoP_Shock_T_Index( infoP_Shock_T_sort(:,4)<=P_thrd,4 );
Cell_Shockdown_id.T4              = infoP_Shock_T_Index( infoP_Shock_T_sort(:,5)<=P_thrd,5 );
Cell_Shockdown_id.T5              = infoP_Shock_T_Index( infoP_Shock_T_sort(:,6)<=P_thrd,6 );

Cell_Shockdown_id.T_all           = {};
Cell_Shockdown_id.T_all{1}        = Cell_Shockdown_id.T; 
Cell_Shockdown_id.T_all{2}        = Cell_Shockdown_id.T1; 
Cell_Shockdown_id.T_all{3}        = Cell_Shockdown_id.T2; 
Cell_Shockdown_id.T_all{4}        = Cell_Shockdown_id.T3; 
Cell_Shockdown_id.T_all{5}        = Cell_Shockdown_id.T4; 
Cell_Shockdown_id.T_all{6}        = Cell_Shockdown_id.T5; 

Cell_Shockdown_id.T_index         = [];

Cell_Shockdown_id.T_index(:,1)    = (infoP_Shock_T(:,1) <=P_thrd);
Cell_Shockdown_id.T_index(:,2)    = (infoP_Shock_T(:,2) <=P_thrd);
Cell_Shockdown_id.T_index(:,3)    = (infoP_Shock_T(:,3) <=P_thrd);
Cell_Shockdown_id.T_index(:,4)    = (infoP_Shock_T(:,4) <=P_thrd);
Cell_Shockdown_id.T_index(:,5)    = (infoP_Shock_T(:,5) <=P_thrd);
Cell_Shockdown_id.T_index(:,6)    = (infoP_Shock_T(:,6) <=P_thrd);


for trialidx=1:6

    Cell_Shockdown_id.T_Num(trialidx)        = length(Cell_Shockdown_id.T_all{trialidx}) ;
    trialidx

end


for trialidx                      = 2:5

    Cell_Shockdown_id.T_index(:,7)  = sum(Cell_Shockdown_id.T_index(:,trialidx:trialidx+1),2);
    Cell_Shockdown_id.T_index_Num(trialidx) = sum(  Cell_Shockdown_id.T_index(:,7)==2 );
    trialidx

end

find( Cell_Shockdown_id.T_index(:,7)==3 )



%%-------------Save the id of the Shockup regulated cells--------------------

COMET_Correction.Cell_Shockdown_id  = Cell_Shockdown_id;
COMET_Correction.Cell_Shockdown_idx = ( infoP_Shock_T<=0.05) ;

%%-------------Save the id of the Shockup regulated cells--------------------
color                               = zeros(length(COMET_Correction.Cell_Shockdown_id.T) , 3);
%%---------------------Show the Shockup cells from the first to last trial--

for trialidx = 1:size(Cell_Shockdown_id.T_all,2)
    
   color                                 = []
   Cell_Shockdown_id_brainregion         = []

    for k=1:length(Cell_Shockdown_id.T_all{trialidx})
    
        brain_region_id                  = COMET_Correction.CellRegion_id( Cell_Shockdown_id.T_all{trialidx}(k) );
        brain_region_name                = atlas_new.nameLR{brain_region_id};
        [~,unique_area_name_idx]         = ismember(unique_area_name,brain_region_name,'rows');
        color( k , : )                   = unique_area_colormap( unique_area_name_idx>0 , :);
        Cell_Shockdown_id_brainregion(k) =  brain_region_id;
        
    end
    
    COMET_Correction.Cell_Shockdown_id_brainregion{trialidx} = Cell_Shockdown_id_brainregion;
    
    A_listgood               = (Cell_Shockdown_id.T_all{trialidx});

%     A_good_sparse            = A_neuron_sparse(:,A_neuron_good_idx+1);
%     A_good_sparse            = A_good_sparse(:,A_listgood);
%     A_good_2D_norm           = normalize(A_good_sparse,1,'range');
%     A_color                  = (A_good_2D_norm.*(A_good_2D_norm>0.3))*color;
%     A_color_2D               = reshape( A_color , 1944,1944,3 )*1;
    
    color_mask               = [1,1,1];
    
%     pruned_atlas_mask_color  = reshape( atlas_new.pruned_atlas_mask_crop2COMET, 1944*1944 , 1 )*color_mask;
%     pruned_atlas_mask_color  = reshape( pruned_atlas_mask_color , 1944,1944,3 )*1;
%     A_color_2D2              = (rot90(A_color_2D, -2))+pruned_atlas_mask_color;
%     
%     figure(3000+trialidx)
%     imshow(A_color_2D2)%% not yet correct the coordinates of A in atlas!!!
%     
    N_Rcenter_r2             = COMET_Correction.N_Rcenter;
    
    figure(9500+trialidx)
    imshow(pruned_atlas_mask_color)
    hold on;
    for i=1:length(atlas_new.cortical_neuron_id)       
        
        cell_brainregion_idx = intersect( A_listgood , atlas_new.cortical_neuron_id{i} );
        scatter(N_Rcenter_r2(cell_brainregion_idx,1),N_Rcenter_r2(cell_brainregion_idx,2),8,atlas_new.unique_area{2,i} ,'filled')

        Cell_Shockdown_id.Numbrainregion( trialidx, i) = size(cell_brainregion_idx,1);

        axis equal

    end
    hold off;

 
%     % 用于存储图例条目的句柄
%     h = zeros(1, size(atlas_new.unique_area, 2));
%     % 遍历每种颜色，绘制带有对应颜色的圆圈
%     for i = 1:size(atlas_new.unique_area, 2)
%         colorRGB = atlas_new.unique_area{2, i}; % RGB 值
%         h(i) = plot(NaN,NaN,'o','MarkerEdgeColor',colorRGB,'MarkerFaceColor',colorRGB); % 使用 NaN 来创建图例条目
%     end
%     % 添加图例
%     lgd=legend(h, atlas_new.unique_area(1, :));
%     % 去除图例边框和背景
%     set(lgd, 'Box', 'off', 'Color', 'none');
%     lgd.TextColor = 'black';
%     % 隐藏坐标轴
%     axis off;
%     % exportgraphics(gcf, ['Shockup','_P0.95.png'], 'Resolution', 500);

end    


Color_shockdown_trial        = [255 0     0;
                                255 165   0;
                                0   255   0;
                                0   0   255;
                                139 0   255]/255;    

Marker_shockdown_trial       = ['o'; '+'; '*' ;'x' ;'s'];



figure(10500)

    plot(  sum( C_trace_good( Cell_Shockdown_id.T1,: ),1 )/size(Cell_Shockdown_id.T1,1),'Color', Color_shockup_trial(1,:)  );
    hold on;
    plot(  sum( C_trace_good( Cell_Shockdown_id.T2,: ),1 )/size(Cell_Shockdown_id.T2,1),'Color', Color_shockup_trial(2,:)  );
    hold on;
    plot(  sum( C_trace_good( Cell_Shockdown_id.T3,: ),1 )/size(Cell_Shockdown_id.T3,1),'Color', Color_shockup_trial(3,:)  );
    hold on;
    plot(  sum( C_trace_good( Cell_Shockdown_id.T4,: ),1 )/size(Cell_Shockdown_id.T4,1),'Color', Color_shockup_trial(4,:)  );
    hold on;
    plot(  sum( C_trace_good( Cell_Shockdown_id.T5,: ),1 )/size(Cell_Shockdown_id.T5,1),'Color', Color_shockup_trial(5,:)  );
    hold on; 
        h1 = fill(1:size(C_trace_sort_good,2), speed_COMET_smooth_7006_binary*20,'blue','EdgeColor','none');
    set(h1,'facealpha',.3)
hold on;
%     plot(-speed_COMET_smooth  ,'Color', 'black')
    hold on;
    
    for i=1:numel(miniscopecuebeginFrame)
        % 定义时间间隔
        t1 = miniscopeshockbeginFrame(i); % 开始时间
        t2 = miniscopeshockendFrame(i); % 结束时间
        % 添加阴影
        px = [t1, t2, t2, t1]; % 阴影的 x 坐标
        py = [-2, -2, 12, 12]; % 阴影的 y 坐标
        patch(px, py, 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none'); % 添加阴影
    
        %         line([t1 t1], [-2 min(NumY_show, size(C_trace_sort_cuedown,1))+2], 'Color', 'r', 'LineStyle', '--'); % 添加开始时间的竖直线
        %         line([t2 t2], [-2 min(NumY_show, size(C_trace_sort_cuedown,1))+2], 'Color', 'r', 'LineStyle', '--'); % 添加结束时间的竖直线
    end

%%======================================================================================
for i=1:length(atlas_new.cortical_neuron_id) 
    
    Cell_Shockdown_id.NumCellsall_brainregion(i) = size( atlas_new.cortical_neuron_id{1,i},2);
    
end

Cell_Shockdown_id.percent        = Cell_Shockdown_id.Numbrainregion ./ repmat( Cell_Shockdown_id.NumCellsall_brainregion, size(Cell_Shockdown_id.T_all,2), 1);





for trialidx=2:6


   Cell_cueupshockup_trial(    :,trialidx-1) = Cell_cueup_id.T_index(  :,trialidx) & Cell_Shockup_id.T_index(  :,trialidx );
   Cell_cueupshockdown_trial(  :,trialidx-1) = Cell_cueup_id.T_index(  :,trialidx) & Cell_Shockdown_id.T_index(:,trialidx );
   Cell_cuedownshockup_trial(  :,trialidx-1) = Cell_cuedown_id.T_index(:,trialidx) & Cell_Shockup_id.T_index(  :,trialidx );
   Cell_cuedownshockdown_trial(:,trialidx-1) = Cell_cuedown_id.T_index(:,trialidx) & Cell_Shockdown_id.T_index(:,trialidx );
 
   trialidx

  
end

sum(Cell_cueupshockup_trial,1)
sum(Cell_cueupshockdown_trial,1)
sum(Cell_cuedownshockup_trial,1)
sum(Cell_cuedownshockdown_trial,1)

Cell_cueupshockup_trial_num     = sum(double(Cell_cueupshockup_trial),1);
Cell_cueupshockdown_trial_num   = sum(double(Cell_cueupshockdown_trial),1);
Cell_cuedownshockup_trial_num   = sum(double(Cell_cuedownshockup_trial),1);
Cell_cuedownshockdown_trial_num = sum(double(Cell_cuedownshockdown_trial),1);


Cell_cueupshockup_trial_num_PerCueup   = Cell_cueupshockup_trial_num./Cell_cueup_id.T_Num(2:end);
Cell_cueupshockup_trial_num_PerShockup = Cell_cueupshockup_trial_num./Cell_Shockup_id.T_Num(2:end);

Cell_cueupshockdown_trial_num_PerCueup   = Cell_cueupshockdown_trial_num./Cell_cueup_id.T_Num(2:end);
Cell_cueupshockdown_trial_num_PerShockdown = Cell_cueupshockdown_trial_num./Cell_Shockdown_id.T_Num(2:end);


Cell_cuedownshockup_trial_num_PerCuedown    = Cell_cuedownshockup_trial_num./Cell_cuedown_id.T_Num(2:end);
Cell_cuedownshockup_trial_num_PerShockup     = Cell_cuedownshockup_trial_num./Cell_Shockup_id.T_Num(2:end);


Cell_cuedownshockdown_trial_num_PerCuedown   = Cell_cuedownshockdown_trial_num./Cell_cuedown_id.T_Num(2:end);
Cell_cuedownshockdown_trial_num_PerShockdown = Cell_cuedownshockdown_trial_num./Cell_Shockdown_id.T_Num(2:end);

figure(10600)
bar(Cell_cueupshockup_trial_num_PerCueup        )
figure(10601)
bar(Cell_cueupshockup_trial_num_PerShockup      )

figure(10602)
bar(Cell_cueupshockdown_trial_num_PerCueup      )
figure(10603)
bar(Cell_cueupshockdown_trial_num_PerShockdown  )

figure(10604)
bar(Cell_cuedownshockup_trial_num_PerCuedown    )
figure(10605)
bar(Cell_cuedownshockup_trial_num_PerShockup    )

figure(10606)
bar(Cell_cuedownshockdown_trial_num_PerCuedown  )
figure(10607)
bar(Cell_cuedownshockdown_trial_num_PerShockdown)


figure(10700)
bar( Cell_cueup_id.percent' )
figure(10701)
bar( Cell_cuedown_id.percent' )

figure(10702)
bar( Cell_Shockup_id.percent' )
figure(10703)
bar( Cell_Shockdown_id.percent' )




for trialidx=2:6


    trialidx

    A_idx_cueupshockup        = find( Cell_cueupshockup_trial(:,trialidx-1)     );
    A_idx_cueupshockdown      = find( Cell_cueupshockdown_trial(:,trialidx-1)   );
    A_idx_cuedownshockup      = find( Cell_cuedownshockup_trial(:,trialidx-1)   );
    A_idx_cuedownshockdown    = find( Cell_cuedownshockdown_trial(:,trialidx-1) );

    N_Rcenter_r2              = COMET_Correction.N_Rcenter;
    
    figure(10900)
    subplot(2,3,trialidx-1)
    imshow(pruned_atlas_mask_color)
    hold on;

    for i=1:length(atlas_new.cortical_neuron_id)

        cell_brainregion_idx = intersect( A_idx_cueupshockup , atlas_new.cortical_neuron_id{i} );
        scatter(N_Rcenter_r2(cell_brainregion_idx,1),N_Rcenter_r2(cell_brainregion_idx,2),8,atlas_new.unique_area{2,i} ,'filled')
        axis equal

    end

    figure(10901)
    subplot(2,3,trialidx-1)
    imshow(pruned_atlas_mask_color)
    hold on;

    for i=1:length(atlas_new.cortical_neuron_id)

        cell_brainregion_idx = intersect( A_idx_cueupshockdown , atlas_new.cortical_neuron_id{i} );
        scatter(N_Rcenter_r2(cell_brainregion_idx,1),N_Rcenter_r2(cell_brainregion_idx,2),8,atlas_new.unique_area{2,i} ,'filled')
        axis equal

    end

    figure(10902)
    subplot(2,3,trialidx-1)
    imshow(pruned_atlas_mask_color)
    hold on;

    for i=1:length(atlas_new.cortical_neuron_id)

        cell_brainregion_idx = intersect( A_idx_cuedownshockup , atlas_new.cortical_neuron_id{i} );
        scatter(N_Rcenter_r2(cell_brainregion_idx,1),N_Rcenter_r2(cell_brainregion_idx,2),8,atlas_new.unique_area{2,i} ,'filled')
        axis equal

    end

    figure(10903)
    subplot(2,3,trialidx-1)
    imshow(pruned_atlas_mask_color)
    hold on;

    for i=1:length(atlas_new.cortical_neuron_id)

        cell_brainregion_idx = intersect( A_idx_cuedownshockdown , atlas_new.cortical_neuron_id{i} );
        scatter(N_Rcenter_r2(cell_brainregion_idx,1),N_Rcenter_r2(cell_brainregion_idx,2),8,atlas_new.unique_area{2,i} ,'filled')
        axis equal

    end


    figure(11000)
    subplot(4,1,1)
    plot(sum( C_trace_good(Cell_cueupshockup_trial(:,trialidx-1),: ),1 )/Cell_cueupshockup_trial_num(trialidx-1) ,'Color', Color_shockup_trial(trialidx-1,:) )
    hold on; 
         h1 = fill(1:size(C_trace_sort_good,2), speed_COMET_smooth_7006_binary*20,'blue','EdgeColor','none');
    set(h1,'facealpha',.3)
%     plot(speed_COMET_smooth  ,'Color', 'black')
    hold on;
    for i=1:numel(miniscopecuebeginFrame)
        % 定义时间间隔
        t1 = miniscopecuebeginFrame(i); % 开始时间
        t2 = miniscopecueendFrame(i); % 结束时间
        % 添加阴影
        px = [t1, t2, t2, t1]; % 阴影的 x 坐标
        py = [-2, -2, 30, 30]; % 阴影的 y 坐标
        patch(px, py, 'blue', 'FaceAlpha', 0.1, 'EdgeColor', 'none'); % 添加阴影

        hold on
        t1 = miniscopeshockbeginFrame(i); % 开始时间
        t2 = miniscopeshockendFrame(i); % 结束时间
        % 添加阴影
        px = [t1, t2, t2, t1]; % 阴影的 x 坐标
        py = [-2, -2, 30, 30]; % 阴影的 y 坐标
        patch(px, py, 'black', 'FaceAlpha', 0.1, 'EdgeColor', 'none'); % 添加阴影
        %         line([t1 t1], [-2 min(NumY_show, size(C_trace_sort_cuedown,1))+2], 'Color', 'r', 'LineStyle', '--'); % 添加开始时间的竖直线
        %         line([t2 t2], [-2 min(NumY_show, size(C_trace_sort_cuedown,1))+2], 'Color', 'r', 'LineStyle', '--'); % 添加结束时间的竖直线
    end
    hold on

    subplot(4,1,2)
    plot(sum( C_trace_good(Cell_cueupshockdown_trial(:,trialidx-1),: ),1 )/Cell_cueupshockdown_trial_num(trialidx-1),'Color', Color_shockup_trial(trialidx-1,:)  )
    hold on; 
         h2 = fill(1:size(C_trace_sort_good,2), speed_COMET_smooth_7006_binary*20,'blue','EdgeColor','none');
    set(h2,'facealpha',.3)
%     plot(speed_COMET_smooth  ,'Color', 'black')
    hold on;
    for i=1:numel(miniscopecuebeginFrame)
        % 定义时间间隔
        t1 = miniscopecuebeginFrame(i); % 开始时间
        t2 = miniscopecueendFrame(i); % 结束时间
        % 添加阴影
        px = [t1, t2, t2, t1]; % 阴影的 x 坐标
        py = [-2, -2, 30, 30]; % 阴影的 y 坐标
        patch(px, py, 'blue', 'FaceAlpha', 0.1, 'EdgeColor', 'none'); % 添加阴影

        hold on
        t1 = miniscopeshockbeginFrame(i); % 开始时间
        t2 = miniscopeshockendFrame(i); % 结束时间
        % 添加阴影
        px = [t1, t2, t2, t1]; % 阴影的 x 坐标
        py = [-2, -2, 30, 30]; % 阴影的 y 坐标
        patch(px, py, 'black', 'FaceAlpha', 0.1, 'EdgeColor', 'none'); % 添加阴影
        %         line([t1 t1], [-2 min(NumY_show, size(C_trace_sort_cuedown,1))+2], 'Color', 'r', 'LineStyle', '--'); % 添加开始时间的竖直线
        %         line([t2 t2], [-2 min(NumY_show, size(C_trace_sort_cuedown,1))+2], 'Color', 'r', 'LineStyle', '--'); % 添加结束时间的竖直线
    end
    hold on

    subplot(4,1,3)
    plot(sum( C_trace_good(Cell_cuedownshockup_trial(:,trialidx-1),: ),1 )/Cell_cuedownshockup_trial_num(trialidx-1),'Color', Color_shockup_trial(trialidx-1,:)  )
    hold on; 
         h3 = fill(1:size(C_trace_sort_good,2), speed_COMET_smooth_7006_binary*20,'blue','EdgeColor','none');
    set(h3,'facealpha',.3)
%     plot(speed_COMET_smooth  ,'Color', 'black')
    hold on;
    for i=1:numel(miniscopecuebeginFrame)
        % 定义时间间隔
        t1 = miniscopecuebeginFrame(i); % 开始时间
        t2 = miniscopecueendFrame(i); % 结束时间
        % 添加阴影
        px = [t1, t2, t2, t1]; % 阴影的 x 坐标
        py = [-2, -2, 30, 30]; % 阴影的 y 坐标
        patch(px, py, 'blue', 'FaceAlpha', 0.1, 'EdgeColor', 'none'); % 添加阴影

        hold on
        t1 = miniscopeshockbeginFrame(i); % 开始时间
        t2 = miniscopeshockendFrame(i); % 结束时间
        % 添加阴影
        px = [t1, t2, t2, t1]; % 阴影的 x 坐标
        py = [-2, -2, 30, 30]; % 阴影的 y 坐标
        patch(px, py, 'black', 'FaceAlpha', 0.1, 'EdgeColor', 'none'); % 添加阴影
        %         line([t1 t1], [-2 min(NumY_show, size(C_trace_sort_cuedown,1))+2], 'Color', 'r', 'LineStyle', '--'); % 添加开始时间的竖直线
        %         line([t2 t2], [-2 min(NumY_show, size(C_trace_sort_cuedown,1))+2], 'Color', 'r', 'LineStyle', '--'); % 添加结束时间的竖直线
    end
    hold on

    subplot(4,1,4)
    plot(sum( C_trace_good(Cell_cuedownshockdown_trial(:,trialidx-1),: ),1 )/Cell_cuedownshockdown_trial_num(trialidx-1) ,'Color', Color_shockup_trial(trialidx-1,:) )
    hold on; 
         h4 = fill(1:size(C_trace_sort_good,2), speed_COMET_smooth_7006_binary*20,'blue','EdgeColor','none');
    set(h4,'facealpha',.3)
%     plot(speed_COMET_smooth  ,'Color', 'black')
    hold on;
    for i=1:numel(miniscopecuebeginFrame)
        % 定义时间间隔
        t1 = miniscopecuebeginFrame(i); % 开始时间
        t2 = miniscopecueendFrame(i); % 结束时间
        % 添加阴影
        px = [t1, t2, t2, t1]; % 阴影的 x 坐标
        py = [-2, -2, 30, 30]; % 阴影的 y 坐标
        patch(px, py, 'blue', 'FaceAlpha', 0.1, 'EdgeColor', 'none'); % 添加阴影

        hold on
        t1 = miniscopeshockbeginFrame(i); % 开始时间
        t2 = miniscopeshockendFrame(i); % 结束时间
        % 添加阴影
        px = [t1, t2, t2, t1]; % 阴影的 x 坐标
        py = [-2, -2, 30, 30]; % 阴影的 y 坐标
        patch(px, py, 'black', 'FaceAlpha', 0.1, 'EdgeColor', 'none'); % 添加阴影
        %         line([t1 t1], [-2 min(NumY_show, size(C_trace_sort_cuedown,1))+2], 'Color', 'r', 'LineStyle', '--'); % 添加开始时间的竖直线
        %         line([t2 t2], [-2 min(NumY_show, size(C_trace_sort_cuedown,1))+2], 'Color', 'r', 'LineStyle', '--'); % 添加结束时间的竖直线
    end
    hold on



end


%%======================================================================================



COMET_Correction.Cell_cueupshockup_trial     = Cell_cueupshockup_trial;

COMET_Correction.Cell_cueupshockdown_trial   = Cell_cueupshockdown_trial;
 
COMET_Correction.Cell_cuedownshockdown_trial = Cell_cueupshockdown_trial;

COMET_Correction.Cell_cuedownshockup_trial   = Cell_cueupshockdown_trial;




    color     = []   
   

    Cell_cueupShockup_id_brainregion     = {}

    Cell_cueupShockdown_id_brainregion   = {}

    Cell_cuedownShockup_id_brainregion   = {}

    Cell_cuedownShockdown_id_brainregion = {}




for trialidx = 1:5

%-------------------------------------------------------------------------

    brain_region_id                      = COMET_Correction.CellRegion_id( Cell_cueupshockup_trial(:,trialidx) );

    for brainidx = 1: size(brain_region_id,2)

        brain_region_name                = atlas_new.nameLR{brain_region_id(brainidx)};
        [~,unique_area_name_idx]         = ismember(unique_area_name,brain_region_name,'rows');
        color( brainidx , : ,trialidx)   = unique_area_colormap( unique_area_name_idx>0 , :);
        Cell_cueupShockup_id_brainregion.T{trialidx}(brainidx) =  brain_region_id(brainidx);

    end




%-------------------------------------------------------------------------

    brain_region_id                      = COMET_Correction.CellRegion_id( Cell_cueupshockdown_trial(:,trialidx) );

    for brainidx = 1: size(brain_region_id,2)

        brain_region_name                = atlas_new.nameLR{brain_region_id(brainidx)};
        [~,unique_area_name_idx]         = ismember(unique_area_name,brain_region_name,'rows');
        color( brainidx , : ,trialidx)   = unique_area_colormap( unique_area_name_idx>0 , :);
        Cell_cueupShockdown_id_brainregion.T{trialidx}(brainidx) =  brain_region_id(brainidx);

    end

%-------------------------------------------------------------------------

    brain_region_id                      = COMET_Correction.CellRegion_id( Cell_cuedownshockup_trial(:,trialidx) );

    for brainidx = 1: size(brain_region_id,2)

        brain_region_name                = atlas_new.nameLR{brain_region_id(brainidx)};
        [~,unique_area_name_idx]         = ismember(unique_area_name,brain_region_name,'rows');
        color( brainidx , : ,trialidx)   = unique_area_colormap( unique_area_name_idx>0 , :);
        Cell_cuedownShockup_id_brainregion.T{trialidx}(brainidx) =  brain_region_id(brainidx);

    end

%-------------------------------------------------------------------------

    brain_region_id                      = COMET_Correction.CellRegion_id( Cell_cuedownshockdown_trial(:,trialidx) );
  
    for brainidx = 1: size(brain_region_id,2)

        brain_region_name                = atlas_new.nameLR{brain_region_id(brainidx)};
        [~,unique_area_name_idx]         = ismember(unique_area_name,brain_region_name,'rows');
        color( brainidx , :,trialidx )   = unique_area_colormap( unique_area_name_idx>0 , :);
        Cell_cuedownShockdown_id_brainregion.T{trialidx}(brainidx) =  brain_region_id(brainidx);

    end

%-------------------------------------------------------------------------

trialidx


end


unique(COMET_Correction.CellRegion_id)



for brainidx = 1 : max(COMET_Correction.CellRegion_id)

brainidx







    for trialidx=2:6


%         trialidx


        CellRegion_id_temp          = (COMET_Correction.CellRegion_id == brainidx)';

        Cell_cueupshockup_temp      = Cell_cueupshockup_trial(:,trialidx-1); 
        Chosen_cueupshockup_idx     = Cell_cueupshockup_temp & CellRegion_id_temp;
    
        Cell_cueupshockdown_temp    = Cell_cueupshockdown_trial(:,trialidx-1);
        Chosen_cueupshockdown_idx   = Cell_cueupshockdown_temp & CellRegion_id_temp;

        Cell_cuedownshockup_temp    = Cell_cuedownshockup_trial(:,trialidx-1);
        Chosen_cuedownshockup_idx   = Cell_cuedownshockup_temp & CellRegion_id_temp;

        Cell_cuedownshockdown_temp  = Cell_cuedownshockdown_trial(:,trialidx-1);
        Chosen_cuedownshockdown_idx = Cell_cuedownshockdown_temp & CellRegion_id_temp;


        figure(19999)




        f2000 = figure(20000)
        f2000.WindowState = 'maximized';

        Chosen_idx   = []
        subplot(4,1,1)
        plot(sum( C_trace_good(Chosen_cueupshockup_idx,: ),1 )/sum(Chosen_cueupshockup_idx) ,'Color', Color_shockup_trial(trialidx-1,:) )
        for i=1:numel(miniscopecuebeginFrame)
            % 定义时间间隔
            t1 = miniscopecuebeginFrame(i); % 开始时间
            t2 = miniscopecueendFrame(i); % 结束时间
            % 添加阴影
            px = [t1, t2, t2, t1]; % 阴影的 x 坐标
            py = [-2, -2, 30, 30]; % 阴影的 y 坐标
            patch(px, py, 'blue', 'FaceAlpha', 0.1, 'EdgeColor', 'none'); % 添加阴影

            hold on
            t1 = miniscopeshockbeginFrame(i); % 开始时间
            t2 = miniscopeshockendFrame(i); % 结束时间
            % 添加阴影
            px = [t1, t2, t2, t1]; % 阴影的 x 坐标
            py = [-2, -2, 30, 30]; % 阴影的 y 坐标
            patch(px, py, 'black', 'FaceAlpha', 0.1, 'EdgeColor', 'none'); % 添加阴影
            %         line([t1 t1], [-2 min(NumY_show, size(C_trace_sort_cuedown,1))+2], 'Color', 'r', 'LineStyle', '--'); % 添加开始时间的竖直线
            %         line([t2 t2], [-2 min(NumY_show, size(C_trace_sort_cuedown,1))+2], 'Color', 'r', 'LineStyle', '--'); % 添加结束时间的竖直线
        end
        hold on

        subplot(4,1,2)
        plot(sum( C_trace_good(Chosen_cueupshockdown_idx,: ),1 )/sum(Chosen_cueupshockdown_idx),'Color', Color_shockup_trial(trialidx-1,:)  )
        for i=1:numel(miniscopecuebeginFrame)
            % 定义时间间隔
            t1 = miniscopecuebeginFrame(i); % 开始时间
            t2 = miniscopecueendFrame(i); % 结束时间
            % 添加阴影
            px = [t1, t2, t2, t1]; % 阴影的 x 坐标
            py = [-2, -2, 30, 30]; % 阴影的 y 坐标
            patch(px, py, 'blue', 'FaceAlpha', 0.1, 'EdgeColor', 'none'); % 添加阴影

            hold on
            t1 = miniscopeshockbeginFrame(i); % 开始时间
            t2 = miniscopeshockendFrame(i); % 结束时间
            % 添加阴影
            px = [t1, t2, t2, t1]; % 阴影的 x 坐标
            py = [-2, -2, 30, 30]; % 阴影的 y 坐标
            patch(px, py, 'black', 'FaceAlpha', 0.1, 'EdgeColor', 'none'); % 添加阴影
            %         line([t1 t1], [-2 min(NumY_show, size(C_trace_sort_cuedown,1))+2], 'Color', 'r', 'LineStyle', '--'); % 添加开始时间的竖直线
            %         line([t2 t2], [-2 min(NumY_show, size(C_trace_sort_cuedown,1))+2], 'Color', 'r', 'LineStyle', '--'); % 添加结束时间的竖直线
        end
        hold on

        subplot(4,1,3)
        plot(sum( C_trace_good(Chosen_cuedownshockup_idx,: ),1 )/sum(Chosen_cuedownshockup_idx),'Color', Color_shockup_trial(trialidx-1,:)  )
        for i=1:numel(miniscopecuebeginFrame)
            % 定义时间间隔
            t1 = miniscopecuebeginFrame(i); % 开始时间
            t2 = miniscopecueendFrame(i); % 结束时间
            % 添加阴影
            px = [t1, t2, t2, t1]; % 阴影的 x 坐标
            py = [-2, -2, 30, 30]; % 阴影的 y 坐标
            patch(px, py, 'blue', 'FaceAlpha', 0.1, 'EdgeColor', 'none'); % 添加阴影

            hold on
            t1 = miniscopeshockbeginFrame(i); % 开始时间
            t2 = miniscopeshockendFrame(i); % 结束时间
            % 添加阴影
            px = [t1, t2, t2, t1]; % 阴影的 x 坐标
            py = [-2, -2, 30, 30]; % 阴影的 y 坐标
            patch(px, py, 'black', 'FaceAlpha', 0.1, 'EdgeColor', 'none'); % 添加阴影
            %         line([t1 t1], [-2 min(NumY_show, size(C_trace_sort_cuedown,1))+2], 'Color', 'r', 'LineStyle', '--'); % 添加开始时间的竖直线
            %         line([t2 t2], [-2 min(NumY_show, size(C_trace_sort_cuedown,1))+2], 'Color', 'r', 'LineStyle', '--'); % 添加结束时间的竖直线
        end
        hold on

        subplot(4,1,4)
        plot(sum( C_trace_good(Chosen_cuedownshockdown_idx,: ),1 )/sum(Chosen_cuedownshockdown_idx) ,'Color', Color_shockup_trial(trialidx-1,:) )
        for i=1:numel(miniscopecuebeginFrame)
            % 定义时间间隔
            t1 = miniscopecuebeginFrame(i); % 开始时间
            t2 = miniscopecueendFrame(i); % 结束时间
            % 添加阴影
            px = [t1, t2, t2, t1]; % 阴影的 x 坐标
            py = [-2, -2, 30, 30]; % 阴影的 y 坐标
            patch(px, py, 'blue', 'FaceAlpha', 0.1, 'EdgeColor', 'none'); % 添加阴影

            hold on
            t1 = miniscopeshockbeginFrame(i); % 开始时间
            t2 = miniscopeshockendFrame(i); % 结束时间
            % 添加阴影
            px = [t1, t2, t2, t1]; % 阴影的 x 坐标
            py = [-2, -2, 30, 30]; % 阴影的 y 坐标
            patch(px, py, 'black', 'FaceAlpha', 0.1, 'EdgeColor', 'none'); % 添加阴影
            %         line([t1 t1], [-2 min(NumY_show, size(C_trace_sort_cuedown,1))+2], 'Color', 'r', 'LineStyle', '--'); % 添加开始时间的竖直线
            %         line([t2 t2], [-2 min(NumY_show, size(C_trace_sort_cuedown,1))+2], 'Color', 'r', 'LineStyle', '--'); % 添加结束时间的竖直线
        end
        hold on


    end

pause 

hold off

close(f2000)

end


addpath(genpath(path),'-begin')

% Fit with seqNMF
K          = 2;
L          = 50;
lambda     =.05;
 
shg; clf
display('Running seqNMF on simulated data (2 simulated sequences + noise)')

rowrank    = randperm(size(C_trace_unify_good, 1));

[W,H]      = seqNMF(C_trace_unify_good(rowrank(1:500),:),'K',K, 'L', L,'lambda', lambda);
 
%% Look at factors

figure(30000); SimpleWHPlot(W,H); title('SeqNMF reconstruction')
figure(30001); SimpleWHPlot(W,H,X); title('SeqNMF factors, with raw data')





%%======================================================================================



    NumBrainRegion       = size( COMET_Correction.CortexRegion_Activity_raw , 1);
 
    C_trace_Shockup_Avg_Brainregion = []

    for trialidx         = 1:6


        for brainreg     = 1: NumBrainRegion
    

            flg          = find(COMET_Correction.Cell_Shockup_id_brainregion{trialidx}==brainreg);        
            fieldsname   = fields( Cell_Shockup_id);

            temp         = sum( C_trace_good( Cell_Shockup_id.(fieldsname{trialidx})(flg),: ),1) / size(flg,2);

            C_trace_Shockup_NumNeuron_Brainregion(brainreg , trialidx) = size(flg,2);
            C_trace_Shockup_Avg_Brainregion(brainreg , : , trialidx)   = temp;


        end


    end

    sum(C_trace_Shockup_NumNeuron_Brainregion,1)
  

    for brainregion = 1:85


        figure(10600)

        for trialidx=2:6
            plot(  C_trace_Shockup_Avg_Brainregion(brainregion,:,trialidx) , Color_shockup_trial(trialidx-1)  );
            hold on;          
        end  

        for i=1:numel(miniscopecuebeginFrame)
            % 定义时间间隔
            t1 = miniscopeshockbeginFrame(i); % 开始时间
            t2 = miniscopeshockendFrame(i); % 结束时间
            % 添加阴影
            px = [t1, t2, t2, t1]; % 阴影的 x 坐标
            py = [-2, -2, 12, 12]; % 阴影的 y 坐标
            patch(px, py, 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none'); % 添加阴影

            %         line([t1 t1], [-2 min(NumY_show, size(C_trace_sort_cuedown,1))+2], 'Color', 'r', 'LineStyle', '--'); % 添加开始时间的竖直线
            %         line([t2 t2], [-2 min(NumY_show, size(C_trace_sort_cuedown,1))+2], 'Color', 'r', 'LineStyle', '--'); % 添加结束时间的竖直线
        end
       brainregion
%         pause
     hold off


    end

    C_trace_Shockdown_Avg_Brainregion = []

    for trialidx         = 1:6

        for brainreg     = 1: NumBrainRegion
    
            flg          = find(COMET_Correction.Cell_Shockdown_id_brainregion{trialidx}==brainreg);        
            fieldsname   = fields( Cell_Shockdown_id);

            temp         = sum( C_trace_good( Cell_Shockdown_id.(fieldsname{trialidx})(flg),: ),1) / size(flg,2);
            C_trace_Shockdown_NumNeuron_Brainregion(brainreg , trialidx) = size(flg,2);
            C_trace_Shockdown_Avg_Brainregion(brainreg , : , trialidx)   = temp;

        end

    end



   for brainregion      = 1:85


        figure(10700)

        for trialidx=2:6
            plot(  C_trace_Shockdown_Avg_Brainregion(brainregion,:,trialidx) , colorlist(trialidx-1)  );
            hold on;          
        end  

        for i=1:numel(miniscopecuebeginFrame)
            % 定义时间间隔
            t1 = miniscopeshockbeginFrame(i); % 开始时间
            t2 = miniscopeshockendFrame(i); % 结束时间
            % 添加阴影
            px = [t1, t2, t2, t1]; % 阴影的 x 坐标
            py = [-2, -2, 12, 12]; % 阴影的 y 坐标
            patch(px, py, 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none'); % 添加阴影
            %         line([t1 t1], [-2 min(NumY_show, size(C_trace_sort_cuedown,1))+2], 'Color', 'r', 'LineStyle', '--'); % 添加开始时间的竖直线
            %         line([t2 t2], [-2 min(NumY_show, size(C_trace_sort_cuedown,1))+2], 'Color', 'r', 'LineStyle', '--'); % 添加结束时间的竖直线
        end

        hold off

%         CellRegion_id
%         pause

   end



COMET_Correction.Num_neuron =  [];

for brainidx = 1: size(atlas_new.nameLR,2)
 
    COMET_Correction.Num_neuron(brainidx) = size(find( COMET_Correction.CellRegion_id == brainidx),2);
    
end

figure(10800)
bar(C_trace_Shockup_NumNeuron_Brainregion)

figure(10801)
bar(C_trace_Shockup_NumNeuron_Brainregion./repmat(COMET_Correction.Num_neuron', [1 6]))

figure(10802)
bar(C_trace_Shockdown_NumNeuron_Brainregion)

figure(10803)
bar(C_trace_Shockdown_NumNeuron_Brainregion./repmat(COMET_Correction.Num_neuron', [1 6]))

figure(11000)
for trialidx=2:6
    bar(C_trace_Shockup_NumNeuron_Brainregion(:,trialidx)/COMET_Correction.Num_neuron');
    trialidx
    %         pause
end
%%=========================================================================
%%=========================================================================
%Please be very careful about this step!!!!!!!!!!!!!!!
%--------------------------------------------------------------------------
%%------------------the structure if each brain region---------------------


for brainid   = 1:size(atlas_new.idLR,2)
   
    COMET_Correction.TraceSUM(brainid,:)    = sum(C_trace_good(find(COMET_Correction.CellRegion_id==brainid),:),1);
    
    COMET_Correction.uniTraceSUM(brainid,:) = sum(C_trace_unify_good(find(COMET_Correction.CellRegion_id==brainid),:),1);
   
    COMET_Correction.TraceAVG(brainid,:)    = COMET_Correction.TraceSUM(brainid,:)/COMET_Correction.Num_neuron(brainid);
   
    COMET_Correction.uniTraceAVG(brainid,:) = COMET_Correction.uniTraceSUM(brainid,:)/COMET_Correction.Num_neuron(brainid);
  
    COMET_Correction.uniTraceAVG_norm(brainid,:) = ( COMET_Correction.uniTraceAVG(brainid,:)-min(COMET_Correction.uniTraceAVG(brainid,:)) )/(max(COMET_Correction.uniTraceAVG(brainid,:))-min(COMET_Correction.uniTraceAVG(brainid,:)));
    brainid

end



figure(12000)

for brainid = 1:size(atlas_new.idLR,2)

    plot( COMET_Correction.TraceAVG(brainid,:)+brainid*15,'color',Mycolor(brainid,:));

    text(7009, COMET_Correction.TraceAVG(brainid,end)+brainid*15,[num2str(brainid) '_' char(atlas_new.nameLR{brainid})],'Color','black');
   
    hold on;

end


for brainid = 1:size(atlas_new.idLR,2)

    f=figure(12001)
%     f.Name = [num2str(brainid) '_' char(Brain_Struct.label{brainid})];
%     suptitle([num2str(brainid) '_' char(Brain_Struct.label{brainid})]);
    subplot(2,1,1)
    plot(COMET_Correction.TraceAVG(brainid,:),'color',Mycolor(brainid,:));
    xlim([1 size(C_trace_good,2)])

    subplot(2,1,2)
    plot(COMET_Correction.uniTraceAVG(brainid,:),'color',Mycolor(brainid,:));
    xlim([1 size(C_trace_good,2)])
%     text(7009, Brain_Struct.uniTraceSUM(brainid,end)+brainid*15,[num2str(brainid) '_' char(Brain_Struct.label{brainid})],'Color','black');
%     axis off
    sgtitle([num2str(brainid) '  ' char(atlas_new.nameLR{brainid})]);

    hold off;
%     pause

end

%%------------------the structure if each brain region---------------------


figure(13000)
for brainid=1:size(atlas_new.idLR,2)

    plot(COMET_Correction.uniTraceAVG_norm(brainid,:)*10  +  brainid,'Color',Mycolor(brainid,:))
    hold on;
%     plot(Cortical_Regions.celltracesum(brainid,:))
    text(7009, COMET_Correction.uniTraceAVG_norm(brainid,end)*10  +  brainid,[num2str(brainid) '_' char(atlas_new.nameLR{brainid})],'Color','black');
    hold on;
    hold on;
    for i  = 1:numel(miniscopecuebeginFrame)
        % 定义时间间隔
        t1 = miniscopecuebeginFrame(i); % 开始时间
        t2 = miniscopecueendFrame(i); % 结束时间
        % 添加阴影
        px = [t1, t2, t2, t1]; % 阴影的 x 坐标
        py = [-2, -2, min(200, size(COMET_Correction.uniTraceAVG_norm(brainid,:),1))+10+size(atlas_new.idLR,2), min(200, size(COMET_Correction.uniTraceAVG_norm(brainid,:),1))+10+size(atlas_new.idLR,2)]; % 阴影的 y 坐标
        patch(px, py, 'r', 'FaceAlpha', 0.01, 'EdgeColor', 'none'); % 添加阴影
        hold on;  
        % 定义时间间隔
        t1 = miniscopeshockbeginFrame(i); % 开始时间
        t2 = miniscopeshockendFrame(i); % 结束时间
        % 添加阴影
        px = [t1, t2, t2, t1]; % 阴影的 x 坐标
        py = [-2, -2, min(200, size(COMET_Correction.uniTraceAVG_norm(brainid,:),1))+10+size(atlas_new.idLR,2), min(200, size(COMET_Correction.uniTraceAVG_norm(brainid,:),1))+10+size(atlas_new.idLR,2)]; % 阴影的 y 坐标
        patch(px, py, 'b', 'FaceAlpha', 0.01, 'EdgeColor', 'none'); % 添加阴影
        hold on;
    end

hold on;

end

% save([path ,'mapped_results.mat'], 'atlas_new','COMET_Correction');


