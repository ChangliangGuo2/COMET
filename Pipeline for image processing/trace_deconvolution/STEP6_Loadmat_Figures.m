%##########################################################################
clear all;
close all;
clc;
Folder             = 'F:\2-Cscope Project-FOV-3mm_NA-0d1\Cscope Exp\Mice_New_LinearTrack\PP367_Ses1_11_09_30\Processed\';
%addpath('','-begin');
%##########################################################################
load([Folder 'Animal_all.mat']);
tempFRb        = [];
ms             = [];
ms.trackLength = animal{1,1}.trackLength;
ms.FR          = animal{1,1}.FR;
ms.numNeurons  = length(animal{1,1}.info);

tempFR         = [ms.FR(:,:,1); ms.FR(end:-1:1,:,2)];
tempFR         = tempFR./repmat( max(tempFR,[],2),1,ms.numNeurons);

tempFRb        = [tempFRb tempFR] ;

[maxloc, locidx] = max(tempFRb, [], 1);%max firing location of each cell, in both left and right trials
[B,Idx]          = sort(locidx);
tempFRb          = tempFRb(:,Idx);
%%
figure(100)
imagesc(tempFRb')