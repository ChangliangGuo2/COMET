path='Z:\trace_fear_conditioning\recall\gzc_rasgrf-ai148d-371\My_V4_Miniscope\';
% load([path,'AMF_despeckle_MC_denoised_8bit_caiman_result.mat'])
load([path,'AMF_despeckle_MC_denoised_8bit_caiman_result_handcured.mat'])
pnr_images=imread([path,'AMF_despeckle_MC_denoised_8bitpnr_images_sigma2.png']);
pnr_images=double(pnr_images);
cn_images=imread([path,'AMF_despeckle_MC_denoised_8bit_correlation_images_sigma2.png']);
cn_images=double(cn_images);
% pnr_cn=pnr_images;
pnr_cn=pnr_images.*cn_images;
block_size = [400, 400];
overlap = [20, 20]; %
bigimage_cure_flage=1;
smallimage_cure_flage=0;
[height, width] = size(pnr_images);

num_blocks = ceil([height/block_size(1), width/block_size(2)]);

current_block = [1, 1];
mag=2.5;
% figure('units', 'normalized', 'position', [0 0 1 1]);
bad_id_list = [];
x_CoM = cellfun(@(s) s.CoM(1), coordinates(A_neuron_good_idx+1));
y_CoM = cellfun(@(s) s.CoM(2), coordinates(A_neuron_good_idx+1));


if bigimage_cure_flage
img2show=pnr_cn;
imshow(img2show, [],'border','tight','initialmagnification','fit')
% set (gcf,'Position',[500,50,mag*size(img2show,2),mag*size(img2show,1)])
axis image;
drawnow;
hold on
for i = 1:numel(coordinates(A_neuron_good_idx+1))
    idx=A_neuron_good_idx(i)+1;


    % 在这里绘制轮廓
    plot(coordinates{idx}.coordinates(:,1)-1+2, coordinates{idx}.coordinates(:,2)-1+2, 'r');
    plot(coordinates{idx}.CoM(2)-1+2, coordinates{idx}.CoM(1)-1+2, '.');% scatter(coordinates{idx}.CoM(2)-cols(1)+2,coordinates{idx}.CoM(1)-rows(1)+2)


end
% 初始化一个空的索引数组
indices = [];
while true
    h = drawrectangle;
    rect = h.Position;


    % 检测键盘按键
    while true
        key = waitforbuttonpress;
        if key == 0
            selection_type = get(gcf, 'SelectionType');
            if ~strcmp(selection_type, 'normal')
                break;
            end
        elseif key == 1
            c = get(gcf, 'CurrentCharacter');
            if c == ' '  % 如果按下的是空格键
                break;  % 结束循环
            end
        end
    end
    xmin = rect(1);
    ymin = rect(2);
    xmax = xmin + rect(3);
    ymax = ymin + rect(4);
    indices = [indices, find(y_CoM-1+2 > xmin & y_CoM-1+2 < xmax & x_CoM-1+2 > ymin & x_CoM-1+2 < ymax)];
    bad_id_list=unique([bad_id_list, A_neuron_good_idx(indices)], 'stable')

    if key == 0  % 如果检测到鼠标点击
        selection_type = get(gcf, 'SelectionType');
        if strcmp(selection_type, 'alt')  % 如果点击的是鼠标右键
            continue
        end
    elseif key == 1
        c = get(gcf, 'CurrentCharacter');
        set(gcf, 'CurrentCharacter', '0'); % 重置当前字符
        if c == ' '  % 如果按下的是空格键
            break;  % 结束循环
        end
    end



end
close all;
bad_id_list=unique(bad_id_list, 'stable');
A_neuron_good_idx=setdiff(A_neuron_good_idx, bad_id_list, 'stable');
A_neuron_bad_idx=[A_neuron_bad_idx,bad_id_list];
A_neuron_bad_idx=unique(A_neuron_bad_idx, 'stable');
end

if smallimage_cure_flage
bad_id_list = [];
x_CoM = cellfun(@(s) s.CoM(1), coordinates(A_neuron_good_idx+1));
y_CoM = cellfun(@(s) s.CoM(2), coordinates(A_neuron_good_idx+1));
while true

    rows = (current_block(1)-1)*block_size(1) - (current_block(1)-1)*overlap(1) + (1:block_size(1));
    cols = (current_block(2)-1)*block_size(2) - (current_block(2)-1)*overlap(2) + (1:block_size(2));


    rows = rows(rows <= height);
    cols = cols(cols <= width);


    img2show=pnr_cn(rows, cols);
    imshow(img2show, [],'border','tight','initialmagnification','fit')
    disp(['当前处于第 ', num2str(current_block(1)), ' 行，第 ', num2str(current_block(2)), ' 列，一共 ', num2str(num_blocks(1)), ' 行，', num2str(num_blocks(2)), ' 列'])
    set (gcf,'Position',[500,50,mag*size(img2show,2),mag*size(img2show,1)])
    axis image;
    drawnow;
    hold on

    for i = 1:numel(coordinates(A_neuron_good_idx+1))
        idx=A_neuron_good_idx(i)+1;

        if coordinates{idx}.CoM(1)>rows(1) && coordinates{idx}.CoM(1)<rows(end) && coordinates{idx}.CoM(2)>cols(1) && coordinates{idx}.CoM(2)<cols(end)
            % 在这里绘制轮廓
            plot(coordinates{idx}.coordinates(:,1)-cols(1)+2, coordinates{idx}.coordinates(:,2)-rows(1)+2, 'r');
            plot(coordinates{idx}.CoM(2)-cols(1)+2, coordinates{idx}.CoM(1)-rows(1)+2, '.');% scatter(coordinates{idx}.CoM(2)-cols(1)+2,coordinates{idx}.CoM(1)-rows(1)+2)
        end

    end


    % 初始化一个空的索引数组
    indices = [];

    while true
        h = drawrectangle;
        rect = h.Position;


        % 检测键盘按键
        while true
            key = waitforbuttonpress;
            if key == 0
                selection_type = get(gcf, 'SelectionType');
                if ~strcmp(selection_type, 'normal')
                    break;
                end
            elseif key == 1
                c = get(gcf, 'CurrentCharacter');
                if c == ' '  % 如果按下的是空格键
                    break;  % 结束循环
                end
            end
        end
        xmin = rect(1);
        ymin = rect(2);
        xmax = xmin + rect(3);
        ymax = ymin + rect(4);
        indices = [indices, find(y_CoM-cols(1)+2 > xmin & y_CoM-cols(1)+2 < xmax & x_CoM-rows(1)+2 > ymin & x_CoM-rows(1)+2 < ymax)];
        bad_id_list=unique([bad_id_list, A_neuron_good_idx(indices)], 'stable')

        if key == 0  % 如果检测到鼠标点击
            selection_type = get(gcf, 'SelectionType');
            if strcmp(selection_type, 'alt')  % 如果点击的是鼠标右键
                continue
            end
        elseif key == 1
            c = get(gcf, 'CurrentCharacter');
            set(gcf, 'CurrentCharacter', '0'); % 重置当前字符
            if c == ' '  % 如果按下的是空格键
                break;  % 结束循环
            end
        end



    end



    % waitfor(gcf, 'CurrentCharacter', char(32));
    % set(gcf, 'CurrentCharacter', '0'); % 重置当前字符
    close all;

    % img2show=pnr_images(rows, cols);


    current_block(2) = current_block(2) + 1;
    if current_block(2) > num_blocks(2)
        current_block(2) = 1;
        current_block(1) = current_block(1) + 1;
    end


    if current_block(1) > num_blocks(1)
        break;
    end
end




bad_id_list=unique(bad_id_list, 'stable');
A_neuron_good_idx=setdiff(A_neuron_good_idx, bad_id_list, 'stable');
A_neuron_bad_idx=[A_neuron_bad_idx,bad_id_list];
A_neuron_bad_idx=unique(A_neuron_bad_idx, 'stable');
end
%看一下所有从good神经元里面挑出来地bad神经元
% figure
% imshow(pnr_images, [],'border','tight','initialmagnification','fit')
% axis image;
% drawnow;
% hold on
% for i = 1:numel(bad_id_list)
%     idx=bad_id_list(i)+1;
%     plot(coordinates{idx}.coordinates(:,1), coordinates{idx}.coordinates(:,2), 'r');
%     plot(coordinates{idx}.CoM(2), coordinates{idx}.CoM(1), '.');% scatter(coordinates{idx}.CoM(2)-cols(1)+2,coordinates{idx}.CoM(1)-rows(1)+2)
% end



figure
imshow(cn_images.*pnr_images,[],'border','tight','initialmagnification','fit')
colormap hot;
axis image;
drawnow;
hold on
for i = 1:numel(coordinates(A_neuron_good_idx+1))
    idx=A_neuron_good_idx(i)+1;
    plot(coordinates{idx}.coordinates(:,1), coordinates{idx}.coordinates(:,2), 'r','LineWidth', 0.5);
    % patch([coordinates{idx}.coordinates(:,1); NaN], [coordinates{idx}.coordinates(:,2); NaN], 'r', 'EdgeColor', 'r', 'LineWidth', 2, 'EdgeAlpha', 0.5);
    % plot(coordinates{idx}.CoM(2), coordinates{idx}.CoM(1), '.');
end
drawnow;
% exportgraphics(gcf, 'H:\CM2scope_experimental_data\running\gzc_rasgrf-ai148d-370\caiman_analysis\result.png', 'Resolution', 500);

save([path,'AMF_despeckle_MC_denoised_8bit_caiman_result_handcured.mat'], 'A_neuron_bad_idx', 'A_neuron_good_idx', 'A_neuron_sparse', 'coordinates', 'C_trace', 'C_raw','detrended_trace','-v6');









