clc, clear
close all

%% this file enable use to generate clean brain region
%  atlas information from COSMOS. some analysis also learned from COSMOS.
%  return a mapping fucntion, through which we can look up the neurons
%  where it is blongs
%  last update: 2/16/2020. YZ

%% load atals data
atlas = importdata('atlas_top_projection.mat');

%% show the projection view
clean_cortex_outline = atlas.clean_cortex_outline;
figure, imshow(atlas.clean_cortex_outline, [])

pixel_size = 4.4 / (242 - 74);
bregma_dis = 0.5; % in mm
RSPd_MO_point = [229, 240];
bregma_pos(1) = RSPd_MO_point(1);
bregma_pos(2) = RSPd_MO_point(2) - round(bregma_dis  / pixel_size);

pl_acad_pos = [229, 85];

mo_rsp_pos = [256, 241];

save('bregma_point.mat', 'bregma_pos')
%% crop the atlas to get only cortical area
cut_point = [150, 72;
             151, 73;
             201, 74;
             254, 74;
             304, 74;
             304, 73];
% break the connection   
for i = 1 : size(cut_point, 1)
    clean_cortex_outline(cut_point(i, 1), cut_point(i, 2)) = false;
end
figure, imshow(clean_cortex_outline, [])
% throw away the left part
bw_con = bwconncomp(clean_cortex_outline, 8);
for i = 1 : length(bw_con.PixelIdxList)
    bw_con_part_size(i) = length(bw_con.PixelIdxList{i});
end
[~, max_ind] = max(bw_con_part_size);

target_out_line = bw_con.PixelIdxList{max_ind};

clean_cortex_outline2 = zeros(size(clean_cortex_outline));
clean_cortex_outline2(target_out_line) = 1;
figure, imshow(clean_cortex_outline2)

%% generate a mask for top projection
% fill hole
cortex_no_holes = imfill(clean_cortex_outline2, 'holes');
figure, imshow(cortex_no_holes)

% remove those in top_projection
clean_top_projection = atlas.top_projection .* cortex_no_holes;
figure, imshow(log(clean_top_projection), [])
save('clean_top_projection.mat', 'clean_top_projection')
%% generate skeleton with smooth boundary
up_sampled_outline = imresize(clean_cortex_outline2, 5);

% smooth a little bit
windowSize = 20;
kernel = ones(windowSize) / windowSize ^ 2;
blurryImage = conv2(single(up_sampled_outline), kernel, 'same');
binaryImage = blurryImage > 0.05; % Rethreshold
figure, imshow(binaryImage)
%% find closed region
% % Invert the image to form black curves on a white background
% binaryImage = ~binaryImage;
% % Get rid of huge background that touches the border
% binaryImage = imclearborder(binaryImage);
% % Count the objects remaining
% [~, numberOfClosedRegions] = bwlabel(binaryImage)
% 
% 
% % compare with number of top image
% numel(unique(atlas.top_projection))
% % figure, plot(boundary), 
%% output transparent skeleton
to_show = nan(size(binaryImage));
to_show(binaryImage > 0.5) = 0;
% to_show = imresize(to_show, size(clean_cortex_outline));
figure, pcolor(to_show),  shading interp
colormap(gray(1))
axis equal
axis off
hold on, scatter(500, 500, 30, 'r')
savefig('cortex_atlas_smooth.fig')
save('cortical_out_line_resize_5.mat', 'to_show')

%% 
figure, imagesc(to_show)
%% generate rotation, shift and calling 
% source_bregma = [1146, 868];
% source_pl_acad = [1812, 260];
% % source_mo_rsp = [1104, 1148];
% 
% tform = fitgeotrans([source_bregma; source_pl_acad],...
%                     [bregma_pos; pl_acad_pos], 'nonreflectivesimilarity');
% 
% 
% % transformation = computeRigidTransformation( [source_bregma; source_pl_acad],...
% %                 [bregma_pos; pl_acad_point; mo_rsp]);
% 
% % transformation = transformation / transformation(9);
% % load test img
% cap_img = single(loadtiff('test_cap_image.tif'));
% % apply transform
% 
% B = imwarp(cap_img(end : -1 : 1, end : -1 : 1),tform).';
% figure, imshow(B, [])
% 
% save('trans_mat', 'tform');
% % test
% [x,y] = transformPointsForward(tform,source_bregma(1), source_bregma(2))
% [x,y] = transformPointsForward(tform,source_pl_acad(1), source_pl_acad(2))
% [x,y] = transformPointsForward(tform,source_mo_rsp(1), source_mo_rsp(2))
% 
% % function calculate_rigid_trans_matrix(vector_1, vector_2)
%     % generate a rigid shift to map vector_1 to vector_2
%     point_1_1 = vector_1(1, :);
%     point_1_2 = vector_1(2, :);
%     point_2_1 = vector_2(1, :);
%     point_2_2 = vector_2(2, :);
%     
%     % shift
%     tx = point_2_1(1) - point_1_1(1);
%     ty = point_2_1(2) - point_2_1(2);
%     shift_mat = [1, 0, 0;
%                  0, 1, 0;
%                  tx, ty, 1];
%     p3_1 = point_1_1;
%     p3_2 = [tx, ty, 1] * shift_mat;
%     % scale
%     
%     
%     
%     % rotation
%     d1 = sqrt(sum((point_1_1 - point_1_2).^2));
%     d2 = sqrt(sum((point_2_1 - point_2_2).^2));
%     d3 = sqrt(sum((point_2_2 - (point_2_1 - point_1_1) - point_1_2).^2));
%     cos_theta = (d3^2 - d1^2 - d2^2) / 2 / d1 / d2;
%     if point_2_2(1) - (point_2_1(1) - point_1_1(1)) >= point_1_2(1) % counter clock
%         sign_ro = 1;
%     else
%         sign_ro = -1;
%     end
%     rotation_mat = [cos_theta, sign_ro * sqrt(1 - cos_theta.^2);
%                     -sign_ro * sqrt(1 - cos_theta.^2), cos_theta]
% end