%==========================================================================
%
%==========================================================================
clear all
close all
clc
Folder_Behavior  = 'F:\2-Cscope Project-FOV-3mm_NA-0d1\Cscope Exp\Mice_New_LinearTrack\PP367_Ses1_11_09_30\BehavCam_0_Analysis\Matlab Codes\Behav_Data_Prod';
Folder_Mini_Prod = 'F:\2-Cscope Project-FOV-3mm_NA-0d1\Cscope Exp\Mice_New_LinearTrack\PP367_Ses1_11_09_30\Processed';
Folder_Mini_Time = 'F:\2-Cscope Project-FOV-3mm_NA-0d1\Cscope Exp\Mice_New_LinearTrack\PP367_Ses1_11_09_30\Miniscope';
addpath('F:\2-Cscope Project-FOV-3mm_NA-0d1\Cscope Exp\Mice_New_LinearTrack\PP367_Ses1_11_09_30\Processed\Matlab Codes\Guo Pipeline\','-begin');
load ([Folder_Behavior '\behav_g1.mat']);
load ([Folder_Mini_Prod '\A_C_S_T_clean.mat']);
%-------------------------------------------------------------------------
if (sum(diff(behav.time)==0) > 0)
    disp(['Matching time stamps. Count = ' num2str(sum(diff(behav.time)==0))])
    iidxx = find(diff(behav.time)==0);
    behav.time(iidxx+1) = behav.time(iidxx+1)+1;
end % if the time stamp is not increasing, just it.
%-------------------------------------------------------------------------
    time     = behav.time(~isnan(behav.position(:,1))); % Position 1 is the poistion of the linear track
    %----Only keep the time stamps when the position is not nan.
    
    position = behav.position(~isnan(behav.position(:,1)),:);
    %-----Only keep the position when the position is not nan.
    
    behav.position = interp1(time,position,behav.time);
    %-----Interprate the positions according to the whole behavioral time.Now
    %the positions have values along the whole recording.
    
    dt = median(diff(behav.time/1000));
    %-----Get the time interval of the behavioral recording.
    %-----Median value is chosen.
    
    tempPos = behav.position(:,1);
    [YY,I,Y0,LB,UB] = hampel(behav.time/1000,behav.position(:,1),5*dt,2);
    %-----Hampel filter to the input vector x to detect and remove outliers.
    
    behav.position(:,1) = YY';
    %-----Apply the filterd postions to the behavioral positions
    
    dx = [0; diff(behav.position(:,1))];
    dy = [0; diff(behav.position(:,2))];
    %-----Position change frame by frame.
    
    behav.speed = sqrt((dx).^2+(dy).^2)/dt;
    %-----Speed of the moving.
    behav.speed = smooth(behav.speed',ceil(1/dt));
    %-----Smooth the speed and speed' is used instead of speed.
    
    behav.dt    = dt;
    %-----Save the time inteveral.
%--------------------------------------------------------------------------
    ms.numSegments = size(A_good_3D,3);
    ms.firing      = s_2d';
    T              = readtable([Folder_Mini_Time '\timeStamps.csv'],'Range','B2');
    T_spk          = table2array(T); %% timestamp from neural recording
    T_spk          = T_spk(:,1);
    ms.time        = double(T_spk);

%------------The infomation of the neurons
%--------------------------------------------------------------------------
    
%---Handling NaNs
    
    speedThresh = 10; %---Please double check the speed threshold value----
    binSize     = 2;

    figure(1000)
    plot(behav.time/1000,tempPos)
    hold on
    plot(behav.time/1000,behav.position(:,1),'r')
    hold off
    
   
    ms.pos = interp1(behav.time, behav.position(:,1), (ms.time)); %ms.pos = x position at every mscam frame
    %above 2 lines added by Lucia because ms.pos was undefined up
    %to here, not sure where that was supposed to come from
    idx         = 1;
    animal{idx} = [];
    animal{idx}.name = Folder_Mini_Time;
    animal{idx}.pos  = (ms.pos)';

 
    animal{idx}.firing       = ms.firing;
    animal{idx}.time         = ms.time; % unit is millisecond.    
    animal{idx}.A            = A_good_3D;
    animal{idx}.C            = c_2d;
    animal{idx}.Cn           = C_Raw_clean;
            

%%
        animal = animalLinearSFR(animal,10,0.05,2); 
        %animalLinearSFR( animal, speedThresh, occThresh,binSize )
        animal = animalInfoContent(animal,500);

        animal = animalStability(animal,500);

        
        %%
% % %Save animal file
animalName = 'Animal_all';
save([Folder_Mini_Prod '/' animalName '.mat'],'animal','-v7.3')


