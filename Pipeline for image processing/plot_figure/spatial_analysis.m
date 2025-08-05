numOfallNaNs = cellfun(@(s) sum(isnan(s.coordinates(:))), coordinates);
sum(numOfallNaNs>4);
numOfgood_id_NaNs = cellfun(@(s) sum(isnan(s.coordinates(:))), coordinates(A_neuron_good_idx+1));
num1=sum(numOfgood_id_NaNs>4);
id1=A_neuron_good_idx(numOfgood_id_NaNs>4)+1;
figure
% 初始化 good id list 和 bad id list
good_id_list = [];
bad_id_list = [];
total_cells = sum(numOfgood_id_NaNs>4); % 总的细胞数量
for i=1:total_cells
    A=reshape(full(A_neuron_sparse(:,id1(i))),1944,1944);
    CoM=round(coordinates{id1(i)}.CoM);

    subplot(2, 1, 1);
    
    imshow(A(CoM(1)-10:CoM(1)+10,CoM(2)-10:CoM(2)+10),[])
    title(['Cell ', num2str(i), ' of ', num2str(total_cells)]); 
    subplot(2, 1, 2);
    plot(C_trace(id1(i),:));
    title('Trace');
    figPos = get(gcf, 'Position');
    btnWidth = 100; % 按钮的宽度
    btnHeight = 40; % 按钮的高度
    btnX1 = figPos(3) * 0.75 - btnWidth / 2; % 按钮的 x 坐标
    btnY1 = figPos(4) * 0.6 - btnHeight / 2; % 按钮的 y 坐标
    btnX2 = figPos(3) * 0.75 - btnWidth / 2+100; % 按钮的 x 坐标
    btnY2 = figPos(4) * 0.6 - btnHeight / 2; % 按钮的 y 坐标

    if exist('hGood', 'var'), delete(hGood); end
    if exist('hBad', 'var'), delete(hBad); end

    % 创建 good 按钮
    hGood=uicontrol('Style', 'pushbutton', 'String', 'Good',...
        'Position', [btnX1 btnY1 btnWidth btnHeight],...
        'Callback', ['good_id_list = [good_id_list, id1(', num2str(i), ')]; uiresume(gcbf);']);

    % 创建 bad 按钮
    hBad=uicontrol('Style', 'pushbutton', 'String', 'Bad',...
        'Position', [btnX2 btnY2 btnWidth btnHeight],...
        'Callback', ['bad_id_list = [bad_id_list, id1(', num2str(i), ')]; uiresume(gcbf);']);

    % 暂停，等待用户点击按钮
    uiwait(gcf);
end

% 显示 good id list 和 bad id list
disp('Good ID List:');
disp(good_id_list);
disp('Bad ID List:');
disp(bad_id_list);


