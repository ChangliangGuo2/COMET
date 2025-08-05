atlas_mask=ccf_top_projection;
pos=round([window_front, window_top, AP_pos_pixel, ML_pos_pixel]);%images在脑图谱中的位置 ([image front坐标,image top坐标，AP长度，ML长度])unit：pixel
[pruned_atlas_mask,atlas_mask,resize_rectangle_pos,smoothed_branches] = resize_smooth_ccf(pos,[1944,1944],atlas_mask);
figure
imshow(atlas_mask,[])
hold on
rectangle('Position',resize_rectangle_pos ,'EdgeColor', 'r', 'LineWidth', 2);

figure
imshow(pruned_atlas_mask,[])
hold on
rectangle('Position',resize_rectangle_pos ,'EdgeColor', 'r', 'LineWidth', 2);

figure
imshow(atlas_mask(round(resize_rectangle_pos(2):resize_rectangle_pos(2)+resize_rectangle_pos(4)),round(resize_rectangle_pos(1):resize_rectangle_pos(1)+resize_rectangle_pos(3))),[])

figure



figure
hold on
for i=1:numel(smoothed_branches)
    scatter(smoothed_branches{i}(:,2),smoothed_branches{i}(:,1),4,'k', 'filled')
end
% rectangle('Position', resize_rectangle_pos, 'EdgeColor', 'r', 'LineWidth', 2, 'LineStyle', '--');
axis off
% set(gcf, 'Color', 'none');

% sum=0;
% for i=1:length(cortical_neuron_id)
%     %     sum=sum+length(cortical_neuron_id{i});
%     scatter(1944-center(cortical_neuron_id{i},2)+resize_rectangle_pos(2)+100,1944-center(cortical_neuron_id{i},1)+resize_rectangle_pos(1)+50,30,unique_area{2,i} ,'filled')
% end
for i=1:length(cortical_neuron_id)
    %     sum=sum+length(cortical_neuron_id{i});
    scatter(scale(1)*(N_center_CCF(cortical_neuron_id{i},2)-pos(2))+resize_rectangle_pos(2),scale(2)*(N_center_CCF(cortical_neuron_id{i},1)-pos(1))+resize_rectangle_pos(1),13,unique_area{2,i} ,'filled')
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
exportgraphics(gcf, 'D:\desktop\idea\8mm LFOV\Atlas\legend_figure.png', 'Resolution', 500);











