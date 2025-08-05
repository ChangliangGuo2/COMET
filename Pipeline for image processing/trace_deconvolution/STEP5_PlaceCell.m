%==========================================================================
%
%
%
%==========================================================================
% function [Conseq] = PlaceCell(animalI,animalS)
%Gets the Info results(animalI) and the Stability results(animalS). Returns the indexes
%of cells above the stability and info criteria, who havee 5 or more consequetive bins with firing rate over
%the 95th percentile (a vector).
% wroks when youjust need to run one session
%##########################################################################
clear all;
close all;
clc;
Folder             = 'F:\2-Cscope Project-FOV-3mm_NA-0d1\Cscope Exp\Mice_New_LinearTrack\PP367_Ses1_11_09_30\Processed\';
%addpath('','-begin');
%##########################################################################
load([Folder 'Animal_all.mat']);
%===============================
animalI                  = animal;
animalS                  = animal;
Index                    = animalI{1}.infoP >= 0.95;
IndexS                   = animalS{1}.stabilityCoefP >=0.95 & animalS{1}.stabilityCoef2P >=0.95;
animalI{1}.numNeurons    = length(animalI{1}.infoP);
per                      = zeros(1,animalI{1}.numNeurons);

for i           = 1:animalI{1}.numNeurons
    if Index(i)== 1 & IndexS(i)== 1 %if exceeds both stability and info criteria
        per(i)  = animalI{1}.FR95Per(i);
    end
    
end

leftoright = animalI{1}.FR(:,:,1); %FR 1 and 2 dimensions are for odd and even trials
rightoleft = animalI{1}.FR(:,:,2);
conseq     = zeros(1,animalI{1}.numNeurons);

% might need to change matrix directions
% leftoright=leftoright.';
% rightoleft=rightoleft.';
for i=1:animalI{1}.numNeurons %for each cell
    countl=0; %how many times the firing rate is over 95 percentile
    countr=0;
    if per(i) ~=0
        abovethr=leftoright(:,i)>= per(i);
        for j=1:length(abovethr)
            if abovethr(j)==1
                countl=countl+1;
                if countl==5
                    break
                end
                
            else
                countl=0;
            end
        end
        
        abovethr=rightoleft(:,i)>= per(i);
        for j=1:length(abovethr)
            if abovethr(j)==1
                countr=countr+1;
                if countr==5
                    break
                end
            else
                countr=0;
            end
        end
    end
    if countl>=5 || countr>=5 %if is place in one direction: good enough. Also, generally the >95 firing bins are consequtive for more than 5 bins, so it's okay although the code doesn't do that
        conseq(i)=1;
    end
    Conseq     = conseq;
end

%% stats
activecells    = size(animal{1,1}.firing,2);
pcpercentage   = length(find(conseq))/length(conseq)
infoall        = mean(animal{1,1}.info)
stabilityall   = 0.5*(mean(animal{1,1}.stabilityCoef,'omitnan')+mean(animal{1,1}.stabilityCoef2,'omitnan'))
infoplace      = mean(animal{1,1}.info(find(conseq)))
stabilityplace = 0.5*(mean(animal{1,1}.stabilityCoef(find(conseq)))+mean(animal{1,1}.stabilityCoef2(find(conseq))))
%output
output         = [max(animal{1,1}.trialNum);activecells;pcpercentage;infoall;infoplace;stabilityall;stabilityplace].';
%save output and conseq
save('conseq','conseq');
save('output','output');
open output

animal{1,1}.conseq         = conseq;
animal{1,1}.Numtrails      = max(animal{1,1}.trialNum);
animal{1,1}.pcpercentage   = pcpercentage;
animal{1,1}.infoall        = infoall;
animal{1,1}.stabilityall   = stabilityall;
animal{1,1}.infoplace      = infoplace;
animal{1,1}.stabilityplace = stabilityplace;


animalName = 'Animal_all';
save([Folder '/' animalName '.mat'],'animal','-v7.3')

FR_all = rightoleft';
figure(100)
    imshow(FR_all)