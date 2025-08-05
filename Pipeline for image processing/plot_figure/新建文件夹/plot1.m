% close all
% uiopen('C:\Users\PKU\Documents\WeChat Files\wxid_2b6atoc9my7s22\FileStorage\File\2023-12\CCF_distribution.fig',1)
obj = get(gca,'children');
atlas_mask=obj(21).CData;
% % figure;
% % imshow(atlas_mask, []);
% % hold on; % 保持当前图像，使得矩形可以被绘制在上面
% % rectangle('Position', [50 50 100 100], 'EdgeColor', 'r', 'LineWidth', 2);
pos=round(obj(20).Position);
% atlas_mask_sub=atlas_mask(pos(2):pos(2)+pos(4),pos(1):pos(1)+pos(3));
atlas_mask_sub=atlas_mask;
scale=[1944/pos(3),1944/pos(4)];
atlas_mask_sub=imresize(atlas_mask_sub,[1944/pos(3)*size(atlas_mask,1), 1944/pos(4)*size(atlas_mask,2)]);
resize_rectangle_pos(2)=(pos(2))*1944/pos(4);resize_rectangle_pos(4)=pos(4)*1944/pos(4);
resize_rectangle_pos(1)=(pos(1))*1944/pos(3);resize_rectangle_pos(3)=pos(3)*1944/pos(3);
figure
imshow(atlas_mask_sub,[])
hold on
rectangle('Position',resize_rectangle_pos ,'EdgeColor', 'r', 'LineWidth', 2);


% 1. 骨架化
skeleton = bwmorph(atlas_mask_sub, 'thin', Inf);
% skeleton = bwmorph(skeleton,'thin'); %这应该可以移除那些给你带来问题的像素堆
% skeleton = bwmorph(atlas_mask_sub,'thin');
figure;
prunedSkeleton = bwmorph(skeleton, 'spur', 1);
imshow(skeleton, []);
% 找到所有的分支点
branchPts = bwmorph(prunedSkeleton, 'branchpoints');
endPts = bwmorph(prunedSkeleton, 'endpoints');
hold on
[yb,xb] = find(endPts);
plot(xb, yb, 'r*', 'LineWidth', 2);
[branchPts_row,branchPts_col] = find(branchPts);
plot(branchPts_col, branchPts_row, 'r*', 'LineWidth', 2);


% 从骨架中移除分支点以分割分支
se = strel('square', 3); % 创建一个3x3的结构元素
dilatedBranchPoints = imdilate(branchPts, se);
skeletonWithoutBranches = prunedSkeleton & ~dilatedBranchPoints;
figure;
imshow(skeletonWithoutBranches, [])
%
% 标记连通组件
[L, num] = bwlabel(skeletonWithoutBranches);
figure;
imshow(L, [])
hold on
plot(branchPts_col, branchPts_row, 'r*', 'LineWidth', 2);
% 获取端点和分支点坐标
[epY, epX] = find(endPts);
[bpY, bpX] = find(branchPts);
count=0;
for i=1:num
    i
    Skeleton_tem=(L==i);
    [Y1,X1]= find(Skeleton_tem);
    
    endPts_temp = bwmorph(Skeleton_tem, 'endpoints');
    [Y, X] = find(endPts_temp);
    if length(Y)==1
        continue
    else
        count=count+1;
        [row1, col1]=findNeighbours(prunedSkeleton-skeletonWithoutBranches,Y(1),X(1));
        [row11, col11]=findNeighbours(branchPts+endPts,row1, col1);
        
        [row2, col2]=findNeighbours(prunedSkeleton-skeletonWithoutBranches,Y(2),X(2));
        [row22, col22]=findNeighbours(branchPts+endPts,row2, col2);
        if length(row1) > 1
            row1 = row1(1);
            col1 = col1(1);
        end
        if length(row2) > 1
            row2 = row2(1);
            col2 = col2(1);
        end
        if length(row11) > 1
            row11 = row11(1);
            col11 = col11(1);
        end
        if length(row22) > 1
            row22 = row22(1);
            col22 = col22(1);
        end
        branches{count}=orderSkeletonCoordinates([[Y1;row1;row2;row11;row22] [X1;col1;col2;col11;col22]]);
        temp_smooth=smoothSkeletonWithLowess(branches{count},30);
        %         temp_smooth=[[branches{count}(1,1);temp_smooth(:,1);branches{count}(end,1)] [branches{count}(1,2);temp_smooth(:,2);branches{count}(end,2)]];
        smoothed_branches{count}=smoothbranch(temp_smooth);
    end
end




figure
% imshow(prunedSkeleton,[])
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




% [boundaries, labelMatrix] = bwboundaries(prunedSkeleton, 8);
% imshow(prunedSkeleton,[])
% hold on
% for i=1:numel(boundaries)
%     scatter(boundaries{i}(:,2),boundaries{i}(:,1))
%     w = waitforbuttonpress;
% end



function smoothedCurve = smoothbranch(skeletonCoords)
% 提取 X 和 Y 坐标
x = skeletonCoords(:, 1);
y = skeletonCoords(:, 2);

% 创建参数化变量
t = 1:length(x);

% 使用样条插值创建平滑曲线
ppX = spline(t, x);
ppY = spline(t, y);

% 创建更密集的参数化变量，以便进行插值
tt = linspace(1, length(x), 20 * length(x));

% 计算插值后的坐标
smoothX = ppval(ppX, tt);
smoothY = ppval(ppY, tt);

% 组合 X 和 Y 坐标
smoothedCurve = [smoothX', smoothY'];
end
function smoothedCurve = smoothSkeletonWithMovingAverage(skeletonCoords, windowSize)
% 提取 X 和 Y 坐标
x = skeletonCoords(:, 1);
y = skeletonCoords(:, 2);

% 对 X 和 Y 坐标应用移动平均
smoothX = movmean(x, 2,"Endpoints","shrink");
smoothY = movmean(y, windowSize,"Endpoints","shrink");

% 组合平滑后的 X 和 Y 坐标
smoothedCurve = [smoothX, smoothY];
end
function smoothedCurve = smoothSkeletonWithMovingAverage2(skeletonCoords, windowSize)
% 提取 X 和 Y 坐标
x = skeletonCoords(:, 1);
y = skeletonCoords(:, 2);

% 对中间数据应用移动平均
smoothX = movmean(x, windowSize, 'Endpoints', 'discard');
smoothY = movmean(y, windowSize, 'Endpoints', 'discard');

% 保持首尾数据点不变
smoothX(1:floor(windowSize/2)) = x(1:floor(windowSize/2));
smoothX(end-floor(windowSize/2)+1:end) = x(end-floor(windowSize/2)+1:end);
smoothY(1:floor(windowSize/2)) = y(1:floor(windowSize/2));
smoothY(end-floor(windowSize/2)+1:end) = y(end-floor(windowSize/2)+1:end);

% 组合平滑后的 X 和 Y 坐标
smoothedCurve = [smoothX, smoothY];
end

function smoothedCurve = smoothSkeletonWithLowess(skeletonCoords, windowSize)
% 提取 X 和 Y 坐标
x = skeletonCoords(:, 1);
y = skeletonCoords(:, 2);

% 对 X 和 Y 坐标应用 LOWESS 平滑
smoothX = smoothdata(x, 'loess', windowSize);
smoothY = smoothdata(y, 'loess', windowSize);

% 确保首尾数据点不变
smoothX(1) = x(1);
smoothX(end) = x(end);
smoothY(1) = y(1);
smoothY(end) = y(end);

% 组合平滑后的 X 和 Y 坐标
smoothedCurve = [smoothX, smoothY];
end




function [rows, cols] = findNeighbours(skeleton, row, col)
[rows, cols] = find(skeleton(max(row-1, 1):min(row+1, size(skeleton, 1)), ...
    max(col-1, 1):min(col+1, size(skeleton, 2))));
rows = rows + max(row-1, 1) - 1;
cols = cols + max(col-1, 1) - 1;
% 移除中心点
center = (rows == row & cols == col);
rows(center) = [];
cols(center) = [];
end

function orderedCoords = orderSkeletonCoordinates(skeletonCoords)
% 构建图
G = graph;
numPoints = size(skeletonCoords, 1);
for i = 1:numPoints
    G = addnode(G, num2str(i));
end
for i = 1:numPoints
    for j = i+1:numPoints
        if isNeighbor(skeletonCoords(i,:), skeletonCoords(j,:))
            G = addedge(G, num2str(i), num2str(j));
        end
    end
end

% 找到端点
deg = degree(G);
endpointIndices = find(deg == 1);
if length(endpointIndices) ~= 2
    error('More than two endpoints found. The skeleton must have exactly two endpoints.');
end

% 图的遍历
path = shortestpath(G, num2str(endpointIndices(1)), num2str(endpointIndices(2)));

% 提取有序坐标
orderedCoords = skeletonCoords(str2double(path), :);
end

function isNeigh = isNeighbor(pt1, pt2)
% 检查两点是否相邻（8邻域）
isNeigh = max(abs(pt1 - pt2)) <= 1;
end








