trace=s_2d;
normalized_trace = zeros(size(trace));
for i = 1:size(trace, 1) % 遍历每一个细胞
    min_val = min(trace(i, :));
    max_val = max(trace(i, :));
    normalized_trace(i, :) = (trace(i, :) - min_val) / (max_val - min_val);
end
% imshow(normalized_trace(A_neuron_good_idx+1,:),'border','tight')
imshow(normalized_trace,'border','tight')
set(gca, 'DataAspectRatio', [1 2 1])

% normalized_C_trace = zeros(size(C_trace));
% for i = 1:size(C_trace, 1) % 遍历每一个细胞
%     min_val = min(C_trace(i, :));
%     max_val = max(C_trace(i, :));
%     normalized_C_trace(i, :) = (C_trace(i, :) - min_val) / (max_val - min_val);
% end
% imshow(normalized_C_trace(A_neuron_good_idx+1,:),'border','tight')
% % [W,H] = nnmf(C_trace(A_neuron_good_idx+1,:),8);
% for i=1:4
% subplot(4,1,i)
% [sortedValues, sortedIndices] = sort(W(:,i), 'descend'); % 对第一列降序排序
% topIndices = sortedIndices(1:1000); % 获取前 1000 个最大值的索引
% imshow(normalized_C_trace(A_neuron_good_idx(topIndices)+1,:),'border','tight')
% set(gca, 'DataAspectRatio', [1 5 1]);
% end