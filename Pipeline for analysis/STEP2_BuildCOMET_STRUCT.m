%%=========================================================================
%% 导入数据

clc;

clear;

close all;

path                           = 'D:\XM\FC_M94\';

load( [path 'mapped_results.mat'] );

%%=========================================================================
%%-----------------brain regions atlas-------------------------------------
%%-----------------Correct the atlas mask according to the images online--

figure(100)
    imshow(atlas_new.pruned_atlas_mask_crop2COMET)

figure(101)
    imshow(uint8(label2rgb(atlas_new.clean_top_projection_LR2COMET)))

labeledImage                   = atlas_new.clean_top_projection_LR2COMET;

angleRad                       = deg2rad(-COMET_Correction.Angle2rot_rightclock);

COMET_cranial_posterior_center = COMET_Correction.COMET_cranial_posterior_center;

Size_image_Row                 = size(labeledImage,1);

%%-----------------------------Read video for brain region activity analysis

Raw_denoised                   = tiffreadVolume(  [ path 'AMF_despeckle_MC_denoised_8bit.tif']  );

Raw_denoised_correct           = Raw_denoised(:,:,1);

Raw_denoised_correct           = rot90(Raw_denoised_correct,-2);

R                              = [cos(angleRad), -sin(angleRad), 0; sin(angleRad), cos(angleRad), 0; 0, 0, 1];% 计算旋转矩阵

T1                             = [1, 0, -COMET_cranial_posterior_center(1); 0, 1, -COMET_cranial_posterior_center(2); 0, 0, 1];% 将旋转中心移动到原点

T2                             = [1, 0, COMET_cranial_posterior_center(1); 0, 1, COMET_cranial_posterior_center(2); 0, 0, 1];% 将旋转中心移回原来位置

T                              = T2 * R * T1;% 创建仿射变换矩阵

[height, width, ~]             = size(Raw_denoised_correct);% 获取图像尺寸

[X, Y]                         = meshgrid(1:width, 1:height);% 创建网格

coords                         = [X(:) Y(:) ones(numel(X), 1)]';% 应用仿射变换

new_coords                     = T * coords;% 对新坐标四舍五入并转换为整数

new_x                          = round(reshape(new_coords(1,:), height, width));

new_y                          = round(reshape(new_coords(2,:), height, width));

RotatedImage                   = ImageRotation_Speedup(Raw_denoised_correct,X,Y, new_x,new_y);

RotatedImage                   = imtranslate(uint8(RotatedImage),[0, -(COMET_cranial_posterior_center(2)-Size_image_Row/2)],'FillValues',0);

figure(200)
imshow( uint8(double(RotatedImage) + atlas_new.pruned_atlas_mask_crop2COMET*255) )


instrument_size                = [1944, 1944];;

instrument_pixel_size          = COMET_Correction.instrument_pixel_size;

CCF_pixel_size                 = atlas_new.CCF_pixel_size;

FOV_size_AP                    = (instrument_size(1) * instrument_pixel_size); %FOV AP size (unit:mm)

FOV_size_ML                    = (instrument_size(2) * instrument_pixel_size); %FOV ML size (unit:mm)

AP_pos_pixel                   = FOV_size_AP / CCF_pixel_size;  %FOV AP size (unit:pixel)

ML_pos_pixel                   = FOV_size_ML / CCF_pixel_size;  %FOV ML size (unit:pixel)

top_projection3                = atlas_new.clean_top_projection_LR;

P_pos_pixel                    = (instrument_size * instrument_pixel_size) / CCF_pixel_size; 

scale                          = [instrument_size(1)/AP_pos_pixel,instrument_size(2)/ML_pos_pixel];

top_projection3_Scale          = imresize(top_projection3,[scale(1)*size(top_projection3,1), scale(2)*size(top_projection3,2)],"nearest");

figure(1004)
imshow( label2rgb(top_projection3_Scale),'Colormap',parula );

top_projection3_Scale_c        = interp2(X, Y, double(labeledImage), new_x, new_y,'nearest');

top_projection3_Scale_c        = imtranslate(uint8(top_projection3_Scale_c),[0, (COMET_cranial_posterior_center(2)-Size_image_Row/2)],'FillValues',0);

figure(1005)
imshow( label2rgb(top_projection3_Scale_c),'Colormap',parula );


labeledImage_c                 = imtranslate(double(labeledImage),[0, (COMET_cranial_posterior_center(2)-Size_image_Row/2)],'FillValues',0);

angleRad                       = deg2rad(COMET_Correction.Angle2rot_rightclock);

R                              = [cos(angleRad), -sin(angleRad), 0; sin(angleRad), cos(angleRad), 0; 0, 0, 1];

T1                             = [1, 0, -COMET_cranial_posterior_center(1); 0, 1, -COMET_cranial_posterior_center(2); 0, 0, 1];% 将旋转中心移动到原点

T2                             = [1, 0, COMET_cranial_posterior_center(1); 0, 1, COMET_cranial_posterior_center(2); 0, 0, 1];% 将旋转中心移回原来位置

T                              = T2 * R * T1;% 创建仿射变换矩阵

[height, width, ~]             = size(Raw_denoised_correct);% 获取图像尺寸

[X, Y]                         = meshgrid(1:width, 1:height);% 创建网格

coords                         = [X(:) Y(:) ones(numel(X), 1)]';

new_coords                     = T * coords;% 应用仿射变换

new_x                          = round(reshape(new_coords(1,:), height, width));% 对新坐标四舍五入并转换为整数

new_y                          = round(reshape(new_coords(2,:), height, width));

labeledImage_c                 = interp2(X, Y, double(labeledImage_c), new_x, new_y,'nearest');

labeledImage_c(isnan(labeledImage_c)) = 0;

figure(1006)
imshow( Raw_denoised_correct+uint8(labeledImage_c*1) );

brainidx                       = atlas_new.idLR;
atlas_new.NumPixel_region_COMET_c = []


COMET_roi                      = imread([path 'AMF_despeckle_MC_denoised_8bit_roi.png']);

labeledImage_c_roi             = labeledImage_c .* ( COMET_roi>0 );

figure(1007)
imshow( label2rgb(labeledImage_c_roi), 'Colormap', parula);


for idx                                          = 1:length(brainidx)

    atlas_new.NumPixel_region_COMET_c(idx)       = sum(sum(labeledImage_c_roi==brainidx(idx)));
    idx

end

figure(1008)
plot(atlas_new.NumPixel_region_COMET,'o')
hold on 
plot(atlas_new.NumPixel_region_COMET_c)

Cortex_region_Activity                           = [];
Cortex_region_uniActivity                        = [];

for frameid                                      = 1:size(Raw_denoised,3)
    for k                                        = 1:length(brainidx) %

        temp                                     = ( labeledImage_c==brainidx(k) );
        Raw_denoised_correct                     = Raw_denoised(:,:,frameid);
        Raw_denoised_correct                     = rot90(Raw_denoised_correct,-2);
      
        temp2                                    = double( Raw_denoised_correct ) .* double( temp );
        Cortex_region_Activity(k,frameid)        = sum( temp2,'all');

        if(atlas_new.NumPixel_region_COMET_c(k))
            Cortex_region_uniActivity(k,frameid) = Cortex_region_Activity(k,frameid)/atlas_new.NumPixel_region_COMET_c(k) ;% %
        else
            Cortex_region_uniActivity(k,frameid) = 0;
        end

    end
    frameid

end
  

COMET_Correction.CortexRegion_Activity_raw       = Cortex_region_Activity;
COMET_Correction.CortexRegion_uniActivity_raw    = Cortex_region_uniActivity;

save([path ,'mapped_results.mat'], 'atlas_new','COMET_Correction');

%%----------
%%----------show the activity of brain regions-----------------------------
figure(105)
for brainid                    = 1:length(brainidx) 
    plot(Cortex_region_Activity(brainid,:));
    hold on;
end

figure(106)
for brainid                    = 1:length(brainidx) 
    plot(Cortex_region_uniActivity(brainid,:));
    hold on;
end

Raw_denoised                   = tiffreadVolume([path 'AMF_despeckle_MC_denoised_8bit_nobk.tif']);
Raw_denoised_correct           = Raw_denoised(:,:,1);
Raw_denoised_correct           = rot90(Raw_denoised_correct,-2);

Cortex_region_Activity_nobk    = [];
Cortex_region_uniActivity_nobk = [];

for frameid                    = 1:size(Raw_denoised,3)
    for k                      = 1:length(brainidx) %

        temp                   = ( labeledImage_c==brainidx(k) );
        Raw_denoised_correct   = Raw_denoised(:,:,frameid);
        Raw_denoised_correct   = rot90(Raw_denoised_correct,-2);
      
        temp2                                         = double( Raw_denoised_correct ) .* double( temp );
        Cortex_region_Activity_nobk(k,frameid)        = sum( temp2,'all');

        if( atlas_new.NumPixel_region_COMET_c(k) )
            Cortex_region_uniActivity_nobk(k,frameid) = Cortex_region_Activity_nobk(k,frameid)/atlas_new.NumPixel_region_COMET_c(k) ;% %
        else
            Cortex_region_uniActivity_nobk(k,frameid) = 0;
        end

    end
    frameid

end
 
COMET_Correction.CortexRegion_Activity_nobk           = Cortex_region_Activity_nobk;
COMET_Correction.CortexRegion_uniActivity_nobk        = Cortex_region_uniActivity_nobk;

save(  [path ,'mapped_results.mat'], 'atlas_new','COMET_Correction'  );

figure(107)
for brainid = 1:length(brainidx) 
    plot(Cortex_region_Activity_nobk(brainid,:));
    hold on;
end

figure(108)
for brainid = 1:length(brainidx) 
    plot(Cortex_region_uniActivity_nobk(brainid,:)+brainid*1.5);
    hold on;
end

% % % 
% % % filename = ['D:\XM\FC_M370\' ,'brain_blood_activity_regions.avi'];
% % % 
% % % v        = VideoWriter(filename);  % v - video writer object
% % % open(v)
% % % 
% % % hf=figure(109)
% % % title(atlas_new.nameLR{brainidx(idx)})
% % % 
% % % for idx = 1:length(brainidx) 
% % % 
% % %     h(1)= subplot(4,1,1)
% % %     plot(Cortex_region_uniActivity(idx,:))
% % %     h(2)= subplot(4,1,2)
% % %     plot(Cortex_region_uniActivity_nobk(idx,:));
% % %     h(3)= subplot(4,1,3)
% % %     plot(COMET_Correction.Cortexregion_unibloodvesselActivity(idx,:));
% % %     subplot(4,1,4)
% % %     imshow(atlas_new.clean_top_projection_LR2COMET==brainidx(idx))
% % %     title(atlas_new.nameLR{brainidx(idx)})
% % %     title(h(1), 'Region uniactivity Raw')
% % %     title(h(2), 'Region uniactivity noBK')
% % %     title(h(3), 'Blood Vessel uniactivity')
% % %     F  = getframe(hf);
% % %     writeVideo(v,F) 
% % % %     imwrite(uint8(F.cdata),['D:\XM\FC_M370\' ,'brain_blood_activit.tif'] ,'WriteMode','append');
% % %     idx
% % % %     pause 
% % % 
% % %     hold off
% % % end
% % % close(v)
%%=========================================================================
%%
