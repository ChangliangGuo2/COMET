pnr_images=imread('H:\CM2scope_experimental_data\running\gzc_rasgrf-ai148d-370\caiman_analysis\AMF_MC_denoised_8bitpnr_images_sigma2.png');
cor_images=imread('H:\CM2scope_experimental_data\running\gzc_rasgrf-ai148d-370\caiman_analysis\AMF_MC_denoised_8bit_correlation_images_sigma2.png');
% image = tiffreadVolume('H:\CM2scope_experimental_data\running\gzc_rasgrf-ai148d-370\caiman_analysis\AMF_MC_denoised_8bit.tif');
idx=A_neuron_good_idx(6090)+1;
CoM=round(coordinates{idx}.CoM);
d1           = designfilt('lowpassfir','PassbandFrequency',0.3,'StopbandFrequency',0.4 );
figure

subplot(2, 2, 1);
imshow(pnr_images(CoM(1)-50:CoM(1)+50,CoM(2)-50:CoM(2)+50),[])
colormap hot;
axis image;
hold on
plot(coordinates{idx}.coordinates(:,1)-coordinates{idx}.CoM(2)+52, coordinates{idx}.coordinates(:,2)- coordinates{idx}.CoM(1)+52, 'r','LineWidth', 0.5);
plot(52,52, '.');

subplot(2, 2, 2);
imshow(cor_images(CoM(1)-50:CoM(1)+50,CoM(2)-50:CoM(2)+50),[])
colormap hot;
axis image;
hold on
plot(coordinates{idx}.coordinates(:,1)-coordinates{idx}.CoM(2)+52, coordinates{idx}.coordinates(:,2)- coordinates{idx}.CoM(1)+52, 'w','LineWidth', 0.5);
plot(52,52, '.');

subplot(2, 2, [3,4]);
plot(filtfilt(d1,C_raw(idx,:)));
hold on
plot(C_trace(idx,:),'r');
title('Trace');
% imshow(pnr_images, [],'border','tight','initialmagnification','fit')
img2show=image(CoM(1)-20+1:CoM(1)+20+1,CoM(2)-20+1:CoM(2)+20+1,:);


figure
% sv = sliceViewer(img2show-min(img2show,[],3));
sv = sliceViewer(img2show);





