% 计算每个脑区的神经元数量
neuron_counts = cellfun(@numel, cortical_neuron_id);
% 按神经元数量排序
[sorted_counts, sort_idx] = sort(neuron_counts, 'ascend'); % 按数量升序排序

% 同时调整脑区名称和颜色信息
sorted_areas = unique_area(1, sort_idx);
sorted_colors = unique_area(2, sort_idx);

% 绘制图表
figure;
bar_handles = zeros(1, length(sorted_counts)); % 存储每个条形图的句柄
for i = 1:length(sorted_counts)
    bar_handles(i) = barh(i, sorted_counts(i), 'FaceColor', sorted_colors{i});
    hold on;
end

% 设置Y轴的标签为脑区名
set(gca, 'ytick', 1:length(sorted_counts), 'yticklabel', sorted_areas);

% 添加图例
% lgd=legend(bar_handles, sorted_areas);
% set(lgd, 'Box', 'off', 'Color', 'none','Location', 'southeast');
% 美化图表
xlabel('Neuron Count', 'FontSize', 15); % 增大X轴标签的字体大小
ylabel('Brain Area', 'FontSize', 15); % 增大Y轴标签的字体大小
title('Neuron Count by Brain Area', 'FontSize', 18);
box off; % 去掉边框
axis tight; % 紧凑轴
grid off; % 关闭网格
exportgraphics(gcf, 'D:\desktop\idea\8mm LFOV\Atlas\数量统计.png', 'Resolution', 500);
