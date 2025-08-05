% create cortical mapping for each neuron
%%=========================================================================
%%
%%
%% 貌似脑图谱配准有几个假设：------------------------------------------------
%% 1.默认脑图谱的pixel size=CCF_pixel_size=4.4 / (242 - 74)=26.2 µm;且不考虑不同小鼠脑的大小差异？
%%---- 默认是一致的，成年雄性小鼠，多只鼠的平均尺寸----------------------------
%% 2.认为在脑图谱中窗后边缘的位置是确定，即每次开窗都保证窗后边缘在脑图谱中的坐标？是的
%%---- 人字骨缝的位置，bregma骨缝，前面大血管，是开颅骨的参考。。。。。。。。。。
%%---- 
%%=========================================================================

clc
clear
close all

path                          = 'D:\XM\FC_M94\';

load( [path 'Matlab codes\color_map.mat'] )

load( [path,'Matlab codes\AMF_despeckle_MC_denoised_8bit_caiman_result_handcured.mat'] )

Cortex_Brgma_ref              = imread( [path 'Matlab codes\Cortex_Brgma_Ref.png'] );


figure(100)
    imshow(uint8(Cortex_Brgma_ref))


%% cortical region estimation
atlas                         = importdata( [path 'Matlab codes\atlas_top_projection_plusN.mat'] );

figure(101)
    imshow( (atlas.clean_cortex_outline) );


%%-------------------------------------------------------------------------


[boundaries, L]               = bwboundaries( atlas.clean_cortex_outline,'holes' );

labeledImage                  = zeros( size(atlas.clean_cortex_outline ) );

mask                          = [];
for k                         = 1:length(boundaries)

    boundary                  = boundaries{k};
  
    mask                      = poly2mask(boundary(:, 2), boundary(:, 1), size(atlas.clean_cortex_outline , 1), size(atlas.clean_cortex_outline , 2));

    labeledImage(mask)        = k;  

end

figure(102)
    imshow( label2rgb(labeledImage),'Colormap',parula  );

figure(103);
    mesh( (labeledImage) );

%%-------------------------------------------------------------------------
%%-------------------------------------------------------------------------

Atlas_PWindow_Edge            = 389; %% this is got from the atlas image

CCF_pixel_size                = 4.4 / (242 - 74);

% CCF pixel size (unit:mm/pixel); don't change %
% how to make sure all mice use the same pixel size, as the brain size may be different in each one.---by Changliang Guo
%%Please give more informaiton about this atlas in this case  

%%=========================================================================

COMET_RefImg                  = 1*single(tiffreadVolume([path,'AMF_despeckle_MC_denoised_8bit.tif'],'PixelRegion', {[1 1 inf], [1 1 inf], [1 1 1]}));

COMET_RefImg                  = rot90(COMET_RefImg,-2);

figure(200)
imshow(  uint8(COMET_RefImg*2) );
imwrite( uint8(COMET_RefImg) ,[path,'AMF_despeckle_MC_denoised_8bit_Ref.tif'] );

%Save the image and find the rotation point in it for later adjustment

COMET_RefImg                   = imread( [[path,'AMF_despeckle_MC_denoised_8bit_Ref.png']]);

%%-----------------------------mark the middle line, anterior and posterior lines in the COMET image

Cranial_dis_inMiniscope        = [182 1795]; % which takes 1795-182+1 = 1614 pixels in miniscope recoding

%%-----------------------------mark the middle line, anterior and posterior lines in the COMET image

addpath('D:\XM\SO_M94','begin')
% angle                        = -1.6; %(unit:°);负值为顺时针旋转 %%How is this angle determined?---Q by Changliang

instrument_pixel_size          = 0.004644; % optical instrument pixel size(unit:mm/pixel);

instrument_size                = [1944, 1944];
Angle2rot_rightclock           = 3.07; % Got by ImageJ. The COMET images need to be rotated right clock by 2.19 degrees.
angleRad                       = deg2rad(-Angle2rot_rightclock);

%-----------------------------935 may need to be changed by mouse
COMET_cranial_posterior_center = [Cranial_dis_inMiniscope(2) 952];% Got by ImageJ('AMF_despeckle_MC_denoised_8bit_Ref.png'). Be carefull about the swap of X and Y

%--------------------------------------------------------------------------
RotatedImage                   = ImageRotation(COMET_RefImg,angleRad,COMET_cranial_posterior_center);
figure(201)
imshow(uint8(RotatedImage*1))
hold on;
xline(COMET_cranial_posterior_center(1),'-',{'posterior cranial edge'},'Color','r')


RotatedImage                   = imtranslate(uint8(RotatedImage),[0, -(COMET_cranial_posterior_center(2)-size(COMET_RefImg,1)/2)],'FillValues',0);
figure(202)
imshow(uint8(RotatedImage*1))
hold on;
xline(COMET_cranial_posterior_center(1),'-',{'posterior cranial edge'},'Color','r')
%--------------------------------------------------------------------------

front_pos                      =  2.8; % cranial window Anterior  edge to the brgma point (unit:mm)
%                                        used to be 3.0 from JianLu, I
%                                        changed it to 2.8 mm by 
  %                                     (Cranial_dis_inMiniscope(2)-Cranial_dis_inMiniscope(1)+1)*instrument_pixel_size
back_pos                       = -4.8; % cranial window Posterior edge to the brgma point (unit:mm)

%%-----------------This is the size of the cranial window in this mouse----

window_size                    = front_pos - back_pos; % window AP size (unit:mm)

%%-----------------This is the field of view of COMET----------------------

FOV_size_AP                    = (instrument_size(1) * instrument_pixel_size); %FOV AP size (unit:mm)
FOV_size_ML                    = (instrument_size(2) * instrument_pixel_size); %FOV ML size (unit:mm)

%%-----------------This is the field of view of COMET----------------------
%%-----------------This is the number of pixels the FOV in Atlas-----------

AP_pos_pixel                   = FOV_size_AP / CCF_pixel_size;  %FOV AP size (unit:pixel)
ML_pos_pixel                   = FOV_size_ML / CCF_pixel_size;  %FOV ML size (unit:pixel)

%%-----------------This is the number of pixels the FOV in Atlas-----------

Virtual_Bregma_AP              = Atlas_PWindow_Edge + round(back_pos/CCF_pixel_size) + 1;
RSPd_MO_point                  = [227, Virtual_Bregma_AP]; %% [纵向，横向] %% Where is this point got?---Q by Changliang Guo

%%---------微调的时候使用上诉参数，可以利用提取到的细胞，上下对称来估计，因为血管有可能歪
%----------我计划通过上面给的2.6 mm与 5 mm来定这个坐标，
%--------------------------------------------------------------------------

bregma_point_CCF(1)            = RSPd_MO_point(1); %% RSPd_MO_point is the cross point of RSPD and MO brain region
bregma_point_CCF(2)            = RSPd_MO_point(2);

%-----bregma_point_CCF(2) should be equal to Virtual_Bregma(1)
%which means the bregma location is 500µm away from RSPs_MO_points in AP
%direction
check_pwindow                  = (Atlas_PWindow_Edge-bregma_point_CCF(2)+1)*CCF_pixel_size;

%--------------------------------------------------------------------------

window_front                   = bregma_point_CCF(2) - (front_pos + (Cranial_dis_inMiniscope(1)*instrument_pixel_size)) / CCF_pixel_size;
window_back                    = bregma_point_CCF(2) - (back_pos  - (size(COMET_RefImg,2)-COMET_cranial_posterior_center(1)+1)*instrument_pixel_size) / CCF_pixel_size;

%------bregma_point_CCF: 191 is the column in the atlas matrix, and 223 is the row in the matrix, 
% so it is correct if we check the coordinates, in which x and y are
% swapped.
%--------------------------------------------------------------------------

window_top                                      = bregma_point_CCF(1) - ML_pos_pixel / 2;
window_bottom                                   = bregma_point_CCF(1) + ML_pos_pixel / 2;
% show a CCF FOV, if you are satisfied with the position
COMET_Correction                                = {};
COMET_Correction.Cranial_dis_inMiniscope        = Cranial_dis_inMiniscope;
COMET_Correction.instrument_pixel_size          = instrument_pixel_size;
COMET_Correction.Angle2rot_rightclock           = Angle2rot_rightclock;
COMET_Correction.COMET_cranial_posterior_center = COMET_cranial_posterior_center;
COMET_Correction.ShiftRow                       = -(COMET_cranial_posterior_center(2)-size(COMET_RefImg,1)/2);


figure(203)
imshow( atlas.clean_cortex_outline )
hold on;
xline(Virtual_Bregma_AP,'White')
hold on;
xline(round(window_front),'White')
hold on;
xline(round(window_back),'White')


figure(204)
atlas_rgb(:,:,1)         = atlas.clean_cortex_outline;
atlas_rgb(:,:,2)         = atlas.clean_cortex_outline;
atlas_rgb(:,:,3)         = 0;
imshow( double(atlas_rgb) )
hold on;
xline(Virtual_Bregma_AP,'-',{'referent point'},'Color','white')
hold on;
xline(Atlas_PWindow_Edge,'-',{'posterior cranial edge'},'Color','r')
hold on;
xline(Atlas_PWindow_Edge-round(window_size/CCF_pixel_size),'-',{'anterior cranial edge_correct'},'Color','r')
hold on;
xline(window_front,'-',{'COMET LEFT'},'Color','white')
hold on;
xline(window_back,'-',{'COMENT RIGHT'},'Color','white')
hold on;
FOV_pos                      = round([window_front, window_top, AP_pos_pixel, ML_pos_pixel]);%images在脑图谱中的位置 ([image front坐标,image top坐标，FOV AP长度，FOV ML长度])unit：pixel
yline(window_top,'-',{'COMET TOP'},'Color','white')
hold on;
yline(window_bottom,'-',{'COMET BOTTOM'},'Color','white')
(window_bottom+window_top)/2

figure(205)
imshow((atlas.top_projection))

clean_cortex_outline_fill    = imfill(atlas.clean_cortex_outline ,'holes');             %填充
figure(206)
imshow((clean_cortex_outline_fill))

se = strel('disk', 1); % 'disk' 类型，半径为 5 的结构元素
clean_cortex_outline_fillerodedMask = imerode(clean_cortex_outline_fill, se);
figure(207);
imshow(clean_cortex_outline_fillerodedMask);



atlas.clean_top_projection   = atlas.top_projection .* clean_cortex_outline_fillerodedMask;
figure(208)
imshow(label2rgb(atlas.clean_top_projection),'Colormap',parula)

clean_top_projection         = zeros(size(atlas.clean_top_projection));
atlas_newmask_unique         = unique(atlas.clean_top_projection);
atlas_new                    = {};
index                        = 1;

for maskid = 1:size(atlas_newmask_unique,1) 
 
   if(atlas_newmask_unique(maskid,1))

       clean_top_projection(atlas.clean_top_projection==atlas_newmask_unique(maskid,1)) = index;
       nameidx = find(  strcmp( atlas.ids, num2str(atlas_newmask_unique(maskid)) )  );
       atlas_new.id(index)   = (index);
       atlas_new.name{index} = atlas.acronyms{nameidx};
       maskid
       index = index+1

   end

end


mask_idcorrect                                     = zeros(size(clean_top_projection));
mask_idcorrect(size(clean_top_projection)/2:end,:) = max(atlas_new.id);
id_AdjustValue                                     = max(atlas_new.id);
mask_idcorrect                                     = mask_idcorrect .*(clean_top_projection>0);

figure(208)
imshow( label2rgb(mask_idcorrect),'Colormap',parula );


clean_top_projection_LR                            = clean_top_projection + mask_idcorrect;

figure(209)
imshow( label2rgb(clean_top_projection_LR),'Colormap',parula );

% clean_top_projection_LR(198,179)     = clean_top_projection_LR(198,180);
% clean_top_projection_LR(304,408)     = clean_top_projection_LR(304,409);
% clean_top_projection_LR(228,335:339) = clean_top_projection_LR(228,333);
% clean_top_projection_LR(229,336:340) = clean_top_projection_LR(228,333);
% clean_top_projection_LR(230,341)     = clean_top_projection_LR(228,333);
% clean_top_projection_LR(243,364)     = clean_top_projection_LR(244,364);
% clean_top_projection_LR(244,366)     = clean_top_projection_LR(244,364);
% clean_top_projection_LR(25,350:351)  = clean_top_projection_LR(26,350);
% clean_top_projection_LR(25,350:351)  = clean_top_projection_LR(26,350);

figure(210)
imshow( label2rgb(clean_top_projection_LR),'Colormap',parula );

figure(211)
imshow(label2rgb(clean_top_projection_LR+255*atlas.clean_cortex_outline),'Colormap',parula)

top_projection                      = zeros(size(clean_top_projection_LR));
atlas_newmask_unique                = unique(clean_top_projection);
index                               = 1;

for maskid                          = 1:size(atlas_newmask_unique,1)
 
   if(atlas_newmask_unique(maskid,1))

       top_projection(clean_top_projection_LR==atlas_newmask_unique(maskid,1)) = index;
       nameidx                   = find(  atlas_new.id==( atlas_newmask_unique(maskid) )   );
       atlas_new.nameLR{index}   = [atlas_new.name{nameidx} '_R'];
       atlas_new.idLR(index)     = (index);
       top_projection(clean_top_projection_LR==atlas_newmask_unique(maskid,1)+id_AdjustValue) = index + 1;
       atlas_new.nameLR{index+1} = [atlas_new.name{nameidx} '_L'];
       atlas_new.idLR(index+1)   = (index+1);
       index                     = index+2; %%分成了左右脑区，所以+2

   end

end


%------------------------leave space for outlier neurons

atlas_new.idLR(index)             = (index);
atlas_new.nameLR{index}           = 'Outliers';

%Added by Changliang Guo at 05312024

%------------------------leave space for outlier neurons

%------------------------Start-Check the brain regions---------------------
% % % 
% % % for idx=1:size(atlas_new.idLR,2)
% % % 
% % %     idx
% % %     atlas_new.nameLR{idx}  
% % %     temp                         = [];
% % %     figure(212)
% % %     set(gcf, 'Resize', 'on', 'Position', [100, 100, 800, 600]);
% % %     temp(:,:,1)                  = double(top_projection==idx);
% % %     temp(:,:,2)                  = double(atlas.clean_cortex_outline);
% % %     temp(:,:,3)                  = 0;
% % % 
% % %     imshow( temp);
% % % 
% % %     pause
% % % 
% % % end

%------------------------End-Check the brain regions-----------------------
figure(213)
imshow( label2rgb(top_projection),'Colormap',parula );


atlas_new.clean_top_projection_LR = top_projection;
atlas_new.clean_cortex_outline    = atlas.clean_cortex_outline ;
atlas_new.clean_top_projection    = atlas.clean_top_projection;
atlas_new.PWindow_Edge            = Atlas_PWindow_Edge; %% this is got from the atlas image
atlas_new.CCF_pixel_size          = CCF_pixel_size;


%%=========================================================================


%%  get center
%==========================================================================
%   The coordinates are from CaImAn Pipeline, in which the values are
%   related to the orientation of the images loaded into CaImAn
%   ------------------------------------------by Changliang Guo, May 4,2024
%==========================================================================
n                                = numel(A_neuron_good_idx);
N_center                         = zeros(n, 2);
for i                            = 1:length(A_neuron_good_idx)
    idx                          = A_neuron_good_idx(i)+1;
    % N_center(i, 2) = instrument_size(1)-coordinates{idx}.CoM(1); % 提取x坐标
    N_center(i, 2)               = coordinates{idx}.CoM(1); % 提取x坐标
    N_center(i, 1)               = coordinates{idx}.CoM(2); % 提取y坐标
end
COMET_Correction.coordinates_old = coordinates;
COMET_Correction.A_bad_idx       = A_neuron_bad_idx;
COMET_Correction.A_good_idx      = A_neuron_good_idx;
COMET_Correction.A_sparse        = A_neuron_sparse;
COMET_Correction.C_raw           = C_raw;
COMET_Correction.C_trace         = C_trace;
% COMET_Correction.C_2d            = c_2d;
COMET_Correction.COMET_RefImg    = COMET_RefImg;
%%------show the distributions of A----------------------------------------
figure(300)
plot(N_center(:,1),N_center(:,2),'o')
xlim([1 1944])
ylim([1 1944])
axis equal
%%-------------------------------------------------------------------------
%%----------------------rotate N_center by cortex angle in image
% angle = 0; % rotate angle, lockwise is positive
%%=========================================================================
%%  The way to rotate the coordinates of cells is likely not correct
%----以下：我不认为以图像中心旋转是对的，坐标的旋转应该与上面COMET图像旋转保持一致。
%%----------Comment by Changliang Guo at May 4, 2024-----------------------
theta                    = (0-180)/180*pi;
N_Rcenter_r1(:,1)        = N_center(:,1)*cos(theta) + N_center(:,2)*sin(theta);
N_Rcenter_r1(:,2)        = N_center(:,2)*cos(theta) - N_center(:,1)*sin(theta);
N_Rcenter_r1             = N_Rcenter_r1 + size(COMET_RefImg,1);
figure(301)
imshow( uint8(RotatedImage) )
hold on;
plot( N_Rcenter_r1(:,1),N_Rcenter_r1(:,2),'o','Color','y' )
xlim([1 1944])
ylim([1 1944])
axis equal
%%--------------Above: rotated the points by 180 degrees-------------------
N_Rcenter_r2             = N_Rcenter_r1-COMET_cranial_posterior_center;
theta                    = (-Angle2rot_rightclock)/180*pi;
N_Rcenter_r2(:,1)        = N_Rcenter_r2(:,1)*cos(theta) + N_Rcenter_r2(:,2)*sin(theta);
N_Rcenter_r2(:,2)        = N_Rcenter_r2(:,2)*cos(theta) - N_Rcenter_r2(:,1)*sin(theta);
N_Rcenter_r2             = N_Rcenter_r2 + COMET_cranial_posterior_center;

figure(302)
imshow( uint8(RotatedImage) )
hold on;
plot( N_Rcenter_r2(:,1),N_Rcenter_r2(:,2),'o','Color','y' )
xlim([1 1944])
ylim([1 1944])
axis equal
%%--------------Above: rotated the points by Angle2rot_rightclock degrees--
N_Rcenter_r2(:,2)       = N_Rcenter_r2(:,2)-(COMET_cranial_posterior_center(2)-size(COMET_RefImg,1)/2);
figure(303)
imshow( uint8(RotatedImage) )
hold on;
plot( N_Rcenter_r2(:,1),N_Rcenter_r2(:,2),'o','Color','y' )
xlim([1 1944])
ylim([1 1944])
axis equal

COMET_Correction.N_Rcenter = N_Rcenter_r2;

%%--------------Above: Shift the points by "COMET_cranial_posterior_center(2)-size(COMET_RefImg,1)/2" pixels
%%=========================================================================
COMET_Ref_points = COMET_cranial_posterior_center;
%%=========================================================================
%%-----------------------
[pruned_atlas_mask,atlas_mask,resize_rectangle_pos,smoothed_branches] = resize_smooth_ccf(FOV_pos,instrument_size,atlas_new.clean_cortex_outline,0);
figure(400)
imshow(uint8(atlas_mask*255))

window_top_Scale                  = round(window_top*CCF_pixel_size/instrument_pixel_size);
window_front_Scale                = round(window_front*CCF_pixel_size/instrument_pixel_size);
ML_pos_pixel_Scale                = round(ML_pos_pixel*CCF_pixel_size/instrument_pixel_size)%should be 1944 in this case
pruned_atlas_mask_crop2COMET      = []
pruned_atlas_mask_crop2COMET      = atlas_mask(window_top_Scale:window_top_Scale+ML_pos_pixel_Scale-1,window_front_Scale:window_front_Scale+ML_pos_pixel_Scale-1);


COMET_Atlas = [];
COMET_Atlas(:,:,1)                = double( pruned_atlas_mask_crop2COMET )*255;
COMET_Atlas(:,:,2)                = double( RotatedImage/1 );
COMET_Atlas(:,:,3)                = 0;

figure(401)
imshow(uint8(double(RotatedImage) + pruned_atlas_mask_crop2COMET*255))
hold on;scatter(N_Rcenter_r2(:,1), N_Rcenter_r2(:,2), 10, 'filled', 'Markerfacecolor', 'r')

atlas_new.pruned_atlas_mask_crop2COMET =pruned_atlas_mask_crop2COMET ;
Scale_Factor                      = size( atlas_new.clean_top_projection_LR)*atlas_new.CCF_pixel_size/COMET_Correction.instrument_pixel_size;

clean_top_projection_LR2Scale     = imresize( atlas_new.clean_top_projection_LR, Scale_Factor,"nearest");
clean_top_projection_LR2COMET     = clean_top_projection_LR2Scale(window_top_Scale:window_top_Scale+ML_pos_pixel_Scale-1,window_front_Scale:window_front_Scale+ML_pos_pixel_Scale-1);

atlas_new.clean_top_projection_LR2COMET =clean_top_projection_LR2COMET ;

brainidx                          = atlas_new.idLR;
atlas_new.NumPixel_region_COMET   = []

for idx = 1:length(brainidx)

    atlas_new.NumPixel_region_Scale(idx) = sum(sum(clean_top_projection_LR2Scale==brainidx(idx)));
    atlas_new.NumPixel_region(idx)       = sum(sum( atlas_new.clean_top_projection_LR==brainidx(idx)));
    atlas_new.NumPixel_region_COMET(idx) = sum(sum( clean_top_projection_LR2COMET==brainidx(idx)));
    idx

end
%%=========================================================================
%% change to AP-RL
img_center             = [instrument_size(1)  / 2, instrument_size(2)  / 2]; % this one is fixed
img_center_in_CCF      = [bregma_point_CCF(1) , (window_front + window_back) / 2];
% up data AP RL positions
N_center_CCF           = zeros(size(N_Rcenter_r2));
xlabel_annoy           = [];
valid_id               = [];
  
% for each of neurons
for i                  = 1 : size(N_Rcenter_r2, 1)

    i

    curr_pos           = N_Rcenter_r2(i, :); % in xyz

    %CCF_loc is the pixel indexs of the cells in atlas mask 
    CCF_loc            = [];
    CCF_loc(2)         = (curr_pos(1) - img_center(1)) * instrument_pixel_size / CCF_pixel_size + img_center_in_CCF(2);
    CCF_loc(1)         = (curr_pos(2) - img_center(2)) * instrument_pixel_size / CCF_pixel_size + img_center_in_CCF(1);
    CCF_loc            = round(CCF_loc);
    N_center_CCF(i, :) = CCF_loc;
    %CCF_loc is the pixel indexs of the cells in atlas mask 

    % query a cortical regions
    id            = atlas_new.clean_top_projection_LR(CCF_loc(1), CCF_loc(2));
    abbr_ind      = find(atlas_new.idLR==id); %# Get all capital letters
    if ~isempty(abbr_ind)
        abbr_ind  = abbr_ind;
    else
        abbr_ind  = size(atlas_new.idLR,2);
    end
    valid_id(i)     = abbr_ind;
    xlabel_annoy{i} = atlas_new.nameLR{abbr_ind};
end

addpath([path 'othercolor'],'begin');
Mycolor           = othercolor('Cat_12');

% find unique brain regions
unique_id         = unique(valid_id);
for i = 1 : length(unique_id)
    cortical_neuron_id{i} = find(unique_id(i) == valid_id);
end
% unique regions
unique_area(1,:) = atlas_new.nameLR(unique_id);
% lable colors in unique_area
for i = 1: (length(unique_area))
    if ~isempty(unique_area{1, i})
        unique_area{2,i} = Mycolor(5*i,:);
    end
end
COMET_Correction.N_center_CCF = N_center_CCF;
atlas_new.valid_id            = valid_id;
atlas_new.cortical_neuron_id  = cortical_neuron_id;
atlas_new.unique_area         = unique_area;
atlas_new.xlabel_annoy        = xlabel_annoy;

%%=========================================================================
%% change to AP-RL
[pruned_atlas_mask,atlas_mask,resize_rectangle_pos,smoothed_branches] = resize_smooth_ccf(FOV_pos,instrument_size,atlas_new.clean_cortex_outline,1);
%% The footprint of extracted neurons distributed in an atlas of origin size.
figure(500)
hold on
for i=1:numel(smoothed_branches)
    scatter(smoothed_branches{i}(:,2),instrument_size(1)-smoothed_branches{i}(:,1),4,'k', 'filled')
end
axis off
for i=1:length(cortical_neuron_id)
    scatter(N_Rcenter_r2(cortical_neuron_id{i},1)+resize_rectangle_pos(1),instrument_size(1)-(N_Rcenter_r2(cortical_neuron_id{i},2)+resize_rectangle_pos(2)),15,unique_area{2,i} ,'filled')
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
figure(600);
imshow(RotatedImage, [])
hold on
for i=1:numel(smoothed_branches)
    xx=smoothed_branches{i}(:,2)-(resize_rectangle_pos(1));yy=smoothed_branches{i}(:,1)-(resize_rectangle_pos(2));
    index=xx<0|xx>instrument_size(1)|yy<0|yy>instrument_size(2);
    xx(index)=[];yy(index)=[];
    scatter(xx,yy,3,'y', 'filled')
end
axis off
for i=1:length(cortical_neuron_id)
    scatter(N_Rcenter_r2(cortical_neuron_id{i},1),N_Rcenter_r2(cortical_neuron_id{i},2),8,unique_area{2,i} ,'filled')
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
% exportgraphics(gcf, [path,'CCF_map2.png'], 'Resolution', 500);

COMET_Correction.window_top_Scale   = window_top_Scale;
COMET_Correction.window_front_Scale = window_front_Scale;
COMET_Correction.ML_pos_pixel_Scale = ML_pos_pixel_Scale;

save([path ,'mapped_results.mat'], 'atlas_new','COMET_Correction')%%考虑放到提取结果的结构体里


