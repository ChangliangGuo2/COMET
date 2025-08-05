% d1           = designfilt('lowpassfir','PassbandFrequency',0.001,'StopbandFrequency',0.01 );
path='H:\CM2scope_experimental_data\trace_fear_conditioning\train\gzc_rasgrf-ai148d-93\My_V4_Miniscope\';
load([path,'AMF_MC_denoised_8bit_caiman_result.mat'])
pnr_images=imread([path,'AMF_MC_denoised_8bitpnr_images_sigma2.png']);
pnr_images=double(pnr_images);
d2           = designfilt('lowpassfir','PassbandFrequency',0.05,'StopbandFrequency',0.1 );
bad_id_list=[];
parfor i         = 1:size(A_neuron_good_idx,2)
    idx=A_neuron_good_idx(i)+1;
    % C_Raw_g1      = C_raw(idx+1,:);

    C_Raw_g2      = C_trace(idx,:);
    % C_Raw_g_flt1 = filtfilt(d1,C_Raw_g2);
    C_Raw_g_flt2 = filtfilt(d2,C_Raw_g2);
    C_Raw_g_flt3  = movmin(C_Raw_g_flt2,500);
    C_Raw_g_flt4  = movmax(C_Raw_g_flt2,500);
    if max(C_Raw_g_flt4)/abs(max(C_Raw_g_flt3))<2
        bad_id_list=[bad_id_list,A_neuron_good_idx(i)]
        % figure
        % plot(C_Raw_g2);hold on;
        % plot(C_Raw_g_flt2);
        % plot(C_Raw_g_flt4);
        % plot(C_Raw_g_flt3);
        % 1
    end

end
% [C_trace_unify,PS] = mapminmax(C_trace,0,1);
% for i = 400:500
%     idx=bad_id_list(i)+1;
%     plot(C_trace_unify(idx,:)*3+i,'k');
%     hold on;
% end
bad_id_list=unique(bad_id_list, 'stable');
A_neuron_good_idx=setdiff(A_neuron_good_idx, bad_id_list, 'stable');
A_neuron_bad_idx=[A_neuron_bad_idx,bad_id_list];
A_neuron_bad_idx=unique(A_neuron_bad_idx, 'stable');

% id1=A_neuron_bad_idx;
% figure
% % 初始化 good id list 和 bad id list
% good_id_list = [];
% bad_id_list = [];
% total_cells = length(id1); % 总的细胞数量
% for i=1:total_cells
%     idx=id1(i)+1;
%     A=reshape(full(A_neuron_sparse(:,idx)),1944,1944);
%     CoM=round(coordinates{idx}.CoM);
% 
%     subplot(2, 2, 1)
%     imshow(A(max(CoM(1)-20,1):min(CoM(1)+20,1944),max(CoM(2)-20,0):min(CoM(2)+20,1944)),[])
%     title(['Cell ', num2str(i), ' of ', num2str(total_cells)]);
%     subplot(2, 2, 2);
%     imshow(pnr_images(max(CoM(1)-20,1):min(CoM(1)+20,1944),max(CoM(2)-20,0):min(CoM(2)+20,1944)),[])
%     hold on
%     plot(coordinates{idx}.coordinates(:,1)-coordinates{idx}.CoM(2)+22, coordinates{idx}.coordinates(:,2)- coordinates{idx}.CoM(1)+22, 'r','LineWidth', 0.5);
%     plot(22,22, '.');
%     subplot(2, 2, [3 4]);
%     plot(C_trace(idx,:));
%     title('Trace');
%     figPos = get(gcf, 'Position');
%     btnWidth = 100; % 按钮的宽度
%     btnHeight = 40; % 按钮的高度
%     btnX1 = figPos(3) * 0.75 - btnWidth / 2; % 按钮的 x 坐标
%     btnY1 = figPos(4) * 0.6 - btnHeight / 2; % 按钮的 y 坐标
%     btnX2 = figPos(3) * 0.75 - btnWidth / 2+100; % 按钮的 x 坐标
%     btnY2 = figPos(4) * 0.6 - btnHeight / 2; % 按钮的 y 坐标
% 
%     if exist('hGood', 'var'), delete(hGood); end
%     if exist('hBad', 'var'), delete(hBad); end
% 
%     % 创建 good 按钮
%     hGood=uicontrol('Style', 'pushbutton', 'String', 'Good',...
%         'Position', [btnX1 btnY1 btnWidth btnHeight],...
%         'Callback', ['good_id_list = [good_id_list, id1(', num2str(i), ')]; uiresume(gcbf);']);
% 
%     % 创建 bad 按钮
%     hBad=uicontrol('Style', 'pushbutton', 'String', 'Bad',...
%         'Position', [btnX2 btnY2 btnWidth btnHeight],...
%         'Callback', ['bad_id_list = [bad_id_list, id1(', num2str(i), ')]; uiresume(gcbf);']);
% 
%     % 暂停，等待用户点击按钮
%     uiwait(gcf);
% end
% 
% % 显示 good id list 和 bad id list
% disp('Good ID List:');
% disp(good_id_list);
% disp('Bad ID List:');
% disp(bad_id_list);