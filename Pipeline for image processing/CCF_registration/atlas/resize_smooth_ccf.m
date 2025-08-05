function [prunedSkeleton,atlas_mask,resize_rectangle_pos,smoothed_branches] = resize_smooth_ccf(pos,img_dim,atlas_mask,smooth_flage)
scale=[img_dim(1)/pos(3),img_dim(2)/pos(4)];
atlas_mask=imresize(atlas_mask,[scale(1)*size(atlas_mask,1), scale(2)*size(atlas_mask,2)]);
resize_rectangle_pos(2)=(pos(2))*scale(2);resize_rectangle_pos(4)=pos(4)*scale(2);
resize_rectangle_pos(1)=(pos(1))*scale(1);resize_rectangle_pos(3)=pos(3)*scale(1);
% figure
% imshow(atlas_mask,[])
% hold on
% rectangle('Position',resize_rectangle_pos ,'EdgeColor', 'r', 'LineWidth', 2);


% 1. 骨架化
skeleton = bwmorph(atlas_mask, 'thin', Inf);
% figure;
prunedSkeleton = bwmorph(skeleton, 'spur', 1);
% imshow(prunedSkeleton, []);
% 找到所有的分支点
branchPts = bwmorph(prunedSkeleton, 'branchpoints');
endPts = bwmorph(prunedSkeleton, 'endpoints');
% hold on
% [yb,xb] = find(endPts);
% plot(xb, yb, 'r*', 'LineWidth', 2);
% [branchPts_row,branchPts_col] = find(branchPts);
% plot(branchPts_col, branchPts_row, 'r*', 'LineWidth', 2);

if smooth_flage
% 从骨架中移除分支点以分割分支
se = strel('square', 3); % 创建一个3x3的结构元素
dilatedBranchPoints = imdilate(branchPts, se);
skeletonWithoutBranches = prunedSkeleton & ~dilatedBranchPoints;
% figure;
% imshow(skeletonWithoutBranches, [])
% 标记连通组件
[L, num] = bwlabel(skeletonWithoutBranches);
% figure;
% imshow(L, [])
% hold on
% plot(branchPts_col, branchPts_row, 'r*', 'LineWidth', 2);
% 获取端点和分支点坐标
[epY, epX] = find(endPts);
[bpY, bpX] = find(branchPts);
count=0;
for i=1:num
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
else
    smoothed_branches=[];
end

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


end