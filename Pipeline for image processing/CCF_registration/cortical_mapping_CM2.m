% create cortical mapping for each neuron
clc
clear
close all
addpath('atlas')
load('atlas\color_map.mat')
%%
path='Z:\trace_fear_conditioning\train\gzc_rasgrf-ai148d-370\My_V4_Miniscope\';
front_pos = 2.8;% cranial window Anterior edge to the brgma point (unit:mm)
back_pos = -5;% cranial window Posterior edge to the brgma point (unit:mm)
RSPd_MO_point = [221, 223]; %% [纵向，横向]
angle=-0.6; %(unit:°)负值为顺时针旋转
instrument_pixel_size = 0.004609; % optical instrument pixel size(unit:mm/pixel);
instrument_size=[1944,1944];
CCF_pixel_size = 4.4 / (242 - 74);% CCF pixel size (unit:mm/pixel); don't change
%%
referece_image=single(tiffreadVolume([path,'AMF_despeckle_MC_denoised_8bit.tif'],'PixelRegion', {[1 1 inf], [1 1 inf], [1 1 1]}));
referece_image=rot90(referece_image,2);
referece_image = imrotate(referece_image, angle, 'bilinear', 'crop');
load([path,'AMF_despeckle_MC_denoised_8bit_caiman_result_handcured.mat'])
%%  get center
n = numel(A_neuron_good_idx);
N_center = zeros(n, 2);
for i = 1:length(A_neuron_good_idx)
    idx=A_neuron_good_idx(i)+1;
    % N_center(i, 2) = instrument_size(1)-coordinates{idx}.CoM(1); % 提取x坐标
    N_center(i, 2) = coordinates{idx}.CoM(1); % 提取x坐标
    N_center(i, 1) = coordinates{idx}.CoM(2); % 提取y坐标
end
%% rotate N_center by cortex angle in image
% angle = 0; % rotate angle, lockwise is positive
N_center=N_center-instrument_size/2;
theta = (angle-180)/180*pi;
N_Rcenter(:,1) = N_center(:,1)*cos(theta) + N_center(:,2)*sin(theta);
N_Rcenter(:,2) = N_center(:,2)*cos(theta) - N_center(:,1)*sin(theta);
N_Rcenter=N_Rcenter+instrument_size/2;
%% cortical region estimation
atlas = importdata('atlas\atlas_top_projection_plusN.mat');

FOV_size_AP = (instrument_size(1) * instrument_pixel_size); %FOV AP size (unit:mm)
FOV_size_ML = (instrument_size(2) * instrument_pixel_size); %FOV ML size (unit:mm)
AP_pos_pixel = FOV_size_AP / CCF_pixel_size;  %FOV AP size (unit:pixel)
ML_pos_pixel = FOV_size_ML / CCF_pixel_size;  %FOV ML size (unit:pixel)

window_size = front_pos - back_pos; % window AP size (unit:mm)
window_margin = (window_size - FOV_size_AP) / 2;

bregma_dis = 0.5; % in mm

bregma_point_CCF(1) = RSPd_MO_point(1);
bregma_point_CCF(2) = RSPd_MO_point(2) - round(bregma_dis  / CCF_pixel_size);

window_front = bregma_point_CCF(2) - (front_pos - window_margin) / CCF_pixel_size;
window_back = bregma_point_CCF(2)  - (back_pos + window_margin) / CCF_pixel_size;

window_top = bregma_point_CCF(1) - ML_pos_pixel / 2;
window_bottom = bregma_point_CCF(1) + ML_pos_pixel / 2;
% show a CCF FOV, if you are satisfied with the position
ccf_top_projection = atlas.clean_cortex_outline;

FOV_pos=round([window_front, window_top, AP_pos_pixel, ML_pos_pixel]);%images在脑图谱中的位置 ([image front坐标,image top坐标，FOV AP长度，FOV ML长度])unit：pixel
[pruned_atlas_mask,atlas_mask,resize_rectangle_pos,smoothed_branches] = resize_smooth_ccf(FOV_pos,instrument_size,ccf_top_projection,0);
atlas_mask_sub=atlas_mask(round(resize_rectangle_pos(2):resize_rectangle_pos(2)+resize_rectangle_pos(4)-1),round(resize_rectangle_pos(1):resize_rectangle_pos(1)+resize_rectangle_pos(3)-1));
figure, imshow((referece_image+atlas_mask_sub*max(referece_image(:)/3))./max(referece_image(:))*3)
hold on, scatter(N_Rcenter(:,1), N_Rcenter(:,2), 10, 'filled', 'Markerfacecolor', 'r')

% figure, imshow(ccf_top_projection, [])
% hold on, rectangle('Position', [window_front, window_top, AP_pos_pixel, ML_pos_pixel], 'EdgeColor','r', 'Rotation', angle);
% hold on, scatter(bregma_point_CCF(2), bregma_point_CCF(1), 20, 'filled', 'Markerfacecolor', 'r')

%% change to AP-RL
img_center = [instrument_size(1)  / 2, instrument_size(2)  / 2]; % this one is fixed
img_center_in_CCF = [bregma_point_CCF(1) , (window_front + window_back) / 2];
% up data AP RL positions

N_center_CCF = zeros(size(N_Rcenter));
xlabel_annoy = [];
valid_id = [];

% for each of neurons
for i = 1 : size(N_Rcenter, 1)
    curr_pos = N_Rcenter(i, :); % in xyz

    %
    CCF_loc = [];
    CCF_loc(2) = (curr_pos(1) - img_center(1)) * instrument_pixel_size / CCF_pixel_size + img_center_in_CCF(2);
    CCF_loc(1) = (curr_pos(2) - img_center(2)) * instrument_pixel_size / CCF_pixel_size + img_center_in_CCF(1);
    CCF_loc = round(CCF_loc);
    N_center_CCF(i, :) = CCF_loc;


    % query a cortical regions
    id = atlas.top_projection(CCF_loc(1), CCF_loc(2));
    abbr_ind = find(strcmp(atlas.ids, num2str(id))); %# Get all capital letters
    if ~isempty(abbr_ind)
        abbr_ind = abbr_ind;
    else
        abbr_ind = 1328;
    end
    valid_id(i) = abbr_ind;
    xlabel_annoy{i} = atlas.acronyms{abbr_ind};
end

% find unique brain regions
unique_id = unique(valid_id);

for i = 1 : length(unique_id)
    cortical_neuron_id{i} = find(unique_id(i) == valid_id);
end
% unique regions
unique_area(1,:) = atlas.acronyms(unique_id);
% lable colors in unique_area
for i = 1: (length(unique_area))
    if ~isempty(unique_area{1, i})
        unique_area{2,i} = color_map_a(i,:);
    end
end

%% The footprint of extracted neurons distributed in an atlas of origin size.
[pruned_atlas_mask,atlas_mask,resize_rectangle_pos,smoothed_branches] = resize_smooth_ccf(FOV_pos,instrument_size,ccf_top_projection,1);

figure
hold on
for i=1:numel(smoothed_branches)
    scatter(smoothed_branches{i}(:,2),instrument_size(1)-smoothed_branches{i}(:,1),4,'k', 'filled')
end
axis off
for i=1:length(cortical_neuron_id)
    scatter(N_Rcenter(cortical_neuron_id{i},1)+resize_rectangle_pos(1),instrument_size(1)-(N_Rcenter(cortical_neuron_id{i},2)+resize_rectangle_pos(2)),15,unique_area{2,i} ,'filled')
end
axis equal


% 用于存储图例条目的句柄
h = zeros(1, size(unique_area, 2));
% 遍历每种颜色，绘制带有对应颜色的圆圈
for i = 1:size(unique_area, 2)
colorRGB = unique_area{2, i}; % RGB 值
h(i) = plot(NaN,NaN,'o','MarkerEdgeColor',colorRGB,'MarkerFaceColor',colorRGB); % 使用 NaN 来创建图例条目
end
% 添加图例
lgd=legend(h, unique_area(1, :));
% 去除图例边框和背景
set(lgd, 'Box', 'off', 'Color', 'none');
% 隐藏坐标轴
axis off;
% exportgraphics(gcf, [path,'CCF_map1.png'], 'Resolution', 500);
%%
figure, imshow(referece_image, [])
hold on
for i=1:numel(smoothed_branches)
    xx=smoothed_branches{i}(:,2)-(resize_rectangle_pos(1));yy=smoothed_branches{i}(:,1)-(resize_rectangle_pos(2));
    index=xx<0|xx>instrument_size(1)|yy<0|yy>instrument_size(2);
    xx(index)=[];yy(index)=[];
    scatter(xx,yy,3,'y', 'filled')
end
axis off
for i=1:length(cortical_neuron_id)
    scatter(N_Rcenter(cortical_neuron_id{i},1),N_Rcenter(cortical_neuron_id{i},2),8,unique_area{2,i} ,'filled')
end
axis equal

% 用于存储图例条目的句柄
h = zeros(1, size(unique_area, 2));
% 遍历每种颜色，绘制带有对应颜色的圆圈
for i = 1:size(unique_area, 2)
colorRGB = unique_area{2, i}; % RGB 值
h(i) = plot(NaN,NaN,'o','MarkerEdgeColor',colorRGB,'MarkerFaceColor',colorRGB); % 使用 NaN 来创建图例条目
end
% 添加图例
lgd=legend(h, unique_area(1, :));
% 设置图例字体颜色为白色
set(lgd, 'TextColor', 'white');
% 去除图例边框和背景
set(lgd, 'Box', 'off', 'Color', 'none');
% 隐藏坐标轴
axis off;
set(gcf, 'units','normalized','outerposition',[0 0 1 1]); 
exportgraphics(gcf, [path,'CCF_map2.png'], 'Resolution', 500);
%%
register_parameter.front_pos=front_pos;
register_parameter.back_pos=back_pos;
register_parameter.RSPd_MO_point=RSPd_MO_point;
register_parameter.angle=angle;
register_parameter.instrument_pixel_size=instrument_pixel_size;
register_parameter.instrument_size=instrument_size;
save([path ,'mapped_results.mat'], 'N_center_CCF', 'cortical_neuron_id', 'xlabel_annoy', 'valid_id', 'unique_area','register_parameter')%%考虑放到提取结果的结构体里



