%% transform DAQ data
clear
clc

tdms_filepath='E:\Behavior_data\20231025\Converted'; %%%%%% Change the ID 
savepath='E:\Behavior_data\20231025\Converted';
% mkdir(savepath);
info = dir([tdms_filepath,'\*.tdms']);
for i = 1:length(info)
    tdms_fullname = fullfile(info(i).folder,info(i).name);

    [data_converted,~,~,~]=convertTDMS(true,tdms_fullname);

    DataSet = table(data_converted.Data.MeasuredData(3).Data,data_converted.Data.MeasuredData(4).Data,data_converted.Data.MeasuredData(5).Data,...
      data_converted.Data.MeasuredData(6).Data, data_converted.Data.MeasuredData(7).Data,data_converted.Data.MeasuredData(8).Data,...
      data_converted.Data.MeasuredData(9).Data,data_converted.Data.MeasuredData(10).Data,'VariableNames',{'Lick','CS1','US1','CS2','US2','VelPul','VelDir','CaMarker'});

    buff=split(info(i).name,".");
    savename=fullfile(savepath,[buff{1},'.mat']);
    save(savename,'DataSet')
  
end

%% plot LickRate in time bin of 1s
clear
clc
warning('off','all')
savepath='E:\Behavior_data\302310_hfTFC_11th\20231104\Converted';
mice=dir([savepath,'\*.mat']);

for j = 1:length(mice) % mice number
    fullname = fullfile(mice(j).folder,mice(j).name);
    load(fullname);
    [~,ID,~] = fileparts(mice(j).name);
%     Result.id{j} = ID;
    
    Lickdata = table2array(DataSet);
    
    for i = 1:height(Lickdata)
                
        if Lickdata(i,6) > 3.2 % Channel6, Pulse  编码器产生的pulse，(512 pulses = 1 circle = 12cm) % speed = Pulse/512*12;
            Lickdata(i,6) = 1;
        else
            Lickdata(i,6) = 0;
        end
        
        if Lickdata(i,7) > 3.2 % Channel7, Direction  编码器产生pulse的方向
            Lickdata(i,7) = 1;
        else
            Lickdata(i,7) = 0;
        end
              
     
    end

    Vel_temp = [0;diff(Lickdata(:,6))];
    a = find(Vel_temp == -1);
    Vel_temp(a) = 1; 
    
       
    for i = 1:length(Vel_temp) %判别逆行运动
        if Lickdata(i,7) == 0 && Vel_temp(i) == 1
            Vel_temp(i) = -1;
        end
    end
 
    
    MoVel = Vel_temp; % Lickdata(:,6)
   
    
end
    
    
 


