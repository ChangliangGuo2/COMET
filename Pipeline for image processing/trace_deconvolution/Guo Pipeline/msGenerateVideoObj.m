function ms = msGenerateVideoObj(dirName)
%MSGENERATEMS Summary of this function goes here
%   Detailed explanation goes here

    MAXFRAMESPERFILE    = 1000; %This is set in the miniscope control software   
    % find avi and csv files
    aviFiles            = dir([dirName '\*.avi']);
    csvFiles            = dir([dirName '\*.csv']);
    
    ms.numFiles         = 0;
    ms.numFrames        = 0;
    ms.vidNum           = [];
    ms.frameNum         = [];
    ms.maxFramesPerFile = MAXFRAMESPERFILE;
    
    %find the total number of relevant video files
    ms.numFiles         = length(aviFiles);
  
    %generate a vidObj for each video file. Also calculate total frames
    for i=1:ms.numFiles

        ms.vidObj{i} = VideoReader([dirName '\' num2str(i-1) '.avi']); %%download "Media Player Codec Pack" in your computer
        ms.vidNum    = [ms.vidNum i*ones(1,ms.vidObj{i}.NumberOfFrames)];
        ms.frameNum  = [ms.frameNum 1:ms.vidObj{i}.NumberOfFrames];
        ms.numFrames = ms.numFrames + ms.vidObj{i}.NumberOfFrames;
        i

    end
    ms.height        = ms.vidObj{1}.Height;
    ms.width         = ms.vidObj{1}.Width;
    
    %read timestamp information
    for i=1:length(csvFiles)
        if strcmp(csvFiles(i).name,'timestamp.csv')

            fileID = fopen([dirName '\' csvFiles(i).name],'r');
            dataArray = textscan(fileID, '%f%f%f%f%[^\n\r]', 'Delimiter', '\t', 'EmptyValue' ,NaN,'HeaderLines' ,1, 'ReturnOnError', false);
            frameNum  = dataArray{:, 1};
            TimeStamp = dataArray{:, 2};
            buffer1   = dataArray{:, 3};
            clearvars dataArray;
            fclose(fileID);

            ms.time    = TimeStamp;
            ms.time(1) = 0;
            ms.maxBufferUsed = max(buffer1);       
          
        end

    end    
   
%     %figure out date and time of recording if that information if available
%     %in folder path
    idx = strfind(dirName, '_');
    idx2 = strfind(dirName,'\');
    if (length(idx) >= 4)
        ms.dateNum = datenum(str2double(dirName((idx(end-2)+1):(idx2(end)-1))), ... %year
            str2double(dirName((idx2(end-1)+1):(idx(end-3)-1))), ... %month
            str2double(dirName((idx(end-3)+1):(idx(end-2)-1))), ... %day
            str2double(dirName((idx2(end)+2):(idx(end-1)-1))), ...%hour
            str2double(dirName((idx(end-1)+2):(idx(end)-1))), ...%minute
            str2double(dirName((idx(end)+2):end)));%second
    end
end

