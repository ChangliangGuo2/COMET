% id1=Cells_list_P(:,2)<0.05;
% id2=Cells_list_P(id1,1);
% A_listgood = A_neuron_good_idx(id2)+1;
% A_good_sparse=A_neuron_sparse(:,A_listgood);
% K = size(length(cortical_neuron_id,2);
% temp = prism;
% color = temp(randi(256,K,1),:);
% A_good_2D_norm = normalize(A_good_sparse,1,'range');
% A_color = (A_good_2D_norm.*(A_good_2D_norm>0.3))*color;
% A_color_2D = reshape( A_color , 1944,1944,3 )*1;
% % figure(100)
% % imshow(A_color_2D) 
% 
% hold on
% 
% id1=Cells_list_P(:,2)>=0.05;
% id2=Cells_list_P(id1,1);
% A_listgood = A_neuron_good_idx(id2)+1;
% A_good_sparse=A_neuron_sparse(:,A_listgood);
% K = size(A_good_sparse,2);
% temp = prism;
% color = ones(K,3);
% A_good_2D_norm = normalize(A_good_sparse,1,'range');
% A_color = (A_good_2D_norm.*(A_good_2D_norm>0.3))*color;
% A_color_2D2 = reshape( A_color , 1944,1944,3 )*1;
% A_color_2D2=A_color_2D2+A_color_2D;
% figure(100)
% imshow(A_color_2D2) 
% imwrite(A_color_2D2,'P_0.05_A.png')

atlas_mask=ccf_top_projection;
pos=round([window_front, window_top, AP_pos_pixel, ML_pos_pixel]);%images在脑图谱中的位置 ([image front坐标,image top坐标，AP长度，ML长度])unit：pixel
[pruned_atlas_mask,atlas_mask,resize_rectangle_pos,smoothed_branches] = resize_smooth_ccf(pos,[1944,1944],atlas_mask);

[~,idx]=sort(Cells_list_P(:,1),'ascend');
Cells_list_P=Cells_list_P(idx,:);
id2=Cells_list_P(Cells_list_P(:,2)>=0.95,1);
A_listgood = A_neuron_good_idx(id2)+1;
A_good_sparse=A_neuron_sparse(:,A_listgood);
color=zeros(length(id2),3);
for i=1:length(cortical_neuron_id)
    for k=1:length(cortical_neuron_id{i})
    if ismember(cortical_neuron_id{i}(k), id2)
    color(id2==cortical_neuron_id{i}(k),:)=color_map_a(i,:);
    end
    end
end
A_good_2D_norm = normalize(A_good_sparse,1,'range');
A_color = (A_good_2D_norm.*(A_good_2D_norm>0.3))*color;
A_color_2D = reshape( A_color , 1944,1944,3 )*1;
% figure(100)
% imshow(A_color_2D) 

hold on

id2=Cells_list_P(Cells_list_P(:,2)<0.95,1);
A_listgood = A_neuron_good_idx(id2)+1;
A_good_sparse=A_neuron_sparse(:,A_listgood);
K = size(A_good_sparse,2);
color = ones(K,3)*1;
A_good_2D_norm = normalize(A_good_sparse,1,'range');
A_color = (A_good_2D_norm.*(A_good_2D_norm>0.3))*color;
A_color_2D2 = reshape( A_color , 1944,1944,3 )*0.15;
A_color_2D2=A_color_2D2+A_color_2D;

color = [1,1,0];
xshift=2;yshift=5;
pruned_atlas_mask_sub=pruned_atlas_mask(round(resize_rectangle_pos(2)-xshift:resize_rectangle_pos(2)+resize_rectangle_pos(4)-xshift-1),round(resize_rectangle_pos(1)-yshift:resize_rectangle_pos(1)+resize_rectangle_pos(3)-yshift-1));
pruned_atlas_mask_color = reshape( pruned_atlas_mask_sub , 1944*1944 , 1 )*color;
pruned_atlas_mask_color = reshape( pruned_atlas_mask_color , 1944,1944,3 )*1;
A_color_2D2=(rot90(A_color_2D2, -2))+pruned_atlas_mask_color;
figure(100)
imshow(A_color_2D2)
hold on
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
lgd.TextColor = 'white';
% 隐藏坐标轴
axis off;
% exportgraphics(gcf, 'H:\CM2scope_experimental_data\running\gzc_rasgrf-ai148d-370\caiman_analysis\P_0.05_A2.png', 'Resolution', 500);
% imwrite(A_color_2D2,'P_0.99_A.png')