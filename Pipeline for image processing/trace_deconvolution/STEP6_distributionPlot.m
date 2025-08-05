
%% Animal Batch
%{
%D:\Hesper\codes_folder\PlaceCell\Original Analysis\Place Cell Analysis
%first column = Information content, second column = stability
% D:\Hesper\codes_folder\PlaceCell\batchdata
animalb={};
%%
load('animal.mat');
load('conseq.mat');
animal{1,1}.conseq=conseq;
%%
animalb{1,1}=animal{1,1};
animalb{1,2}=animal{1,1};
animalb{1,3}=animal{1,1};
animalb{1,4}=animal{1,1};
animalb{1,5}=animal{1,1};
animalb{1,6}=animal{1,1};
animalb{1,7}=animal{1,1};
%%
save('D:\Hesper\codes_folder\PlaceCell\batchdata\Control\E.mat','animalb');
%}
%% 1. Cell sequence plot batch
Controlnames={'A5.mat';'E.mat';'P_A_2nd.mat';'T3.mat';'T4.mat'};
Philpotnames={'Ad.mat';'C.mat';'Fa.mat';'Fb.mat';'P_B_2nd.mat';'P_C_2nd.mat';'A6.mat';'A2.mat'};
%%
tempFRb=[];
names=Controlnames;
for n=1:length(names)
    load(char(names(n)));
    ms=[];
    for d=1:length(animalb)
        if length(animalb{1,d}) > 0
            animal{1,1}=animalb{1,d};
            
            ms.trackLength=animal{1,1}.trackLength;
            ms.FR=animal{1,1}.FR;
            ms.numNeurons=length(animal{1,1}.info);
            
            tempFR = [ms.FR(:,:,1); ms.FR(end:-1:1,:,2)];
            tempFR = tempFR./repmat( max(tempFR,[],2),1,ms.numNeurons);
            
            tempFRb=[tempFRb tempFR] ;
        end
    end
end
[maxloc, locidx] = max(tempFRb, [], 1);%max firing location of each cell, in both left and right trials
[B,Idx]          = sort(locidx);
tempFRb          = tempFRb(:,Idx);
%%
imagesc(tempFRb')
%%
%% 2. correct stability score
names=[Controlnames; Philpotnames];
for n=1:length(names)
    load(char(names(n)));
    ms=[];
    for d=1:length(animalb)
        if isempty(animalb{1,d}) == 0
            animalb{1,d}.stabilitym=(atanh(animalb{1,d}.stabilityCoef) + atanh(animalb{1,d}.stabilityCoef2))/2; %stabilitym = can use as the real stability score
            animalb{1,d}.stabilitymean=mean(animalb{1,d}.stabilitym, 'omitnan');
            animalb{1,d}.stabilitymeanPlace= mean( animalb{1,d}.stabilitym(find(animalb{1,d}.conseq)) );
        end
    end
    save(char(names(n)),'animalb')
end
%% output structure
output=zeros(length(names),2,7);
for n=1:length(names)
    load(char(names(n)));
    for d=1:length(animalb)
        if isempty(animalb{1,d}) == 0
            output(n,1,d)=mean(animalb{1,d}.info);
            output(n,2,d)=mean(animalb{1,d}.info(find(animalb{1,d}.conseq)));
%             output(n,1,d)=animalb{1,d}.stabilitymean;
%             output(n,2,d)=animalb{1,d}.stabilitymeanPlace;
        end
    end
end
%%
output1=reshape(output,13,14);
output1(outputl==0)=nan;
%% Anova
% array1=[];
arraystabilityall=array1(:,1:2:14);
arraystabilityplace=array1(:,2:2:14);
%%
arraystabilityall(arraystabilityall==0)=nan;
[p,tbl,stats] = kruskalwallis(arraystabilityall);
multcompare(stats,'CType','lsd')
%%
% arraystabilityplace(arraystabilityplace==0)=nan;
[p,tbl,stats] = kruskalwallis(arraystabilityplace);
multcompare(stats,'CType','lsd')
%% 3. save stability and information for each cell, separately for control and pilo animal, separately for place and non place cells
%%
ControlAll=[];
PhilAll=[];
ControlPlace=[];
ControlNonplace=[];
PhilPlace=[];
PhilNonplace=[];
%% Control
names=Controlnames;
for n=1:length(names)
    load(char(names(n)));
    for d=1:length(animalb)
        if isempty(animalb{1,d}) == 0
            conseq=animalb{1,d}.conseq;
            CA=[animalb{1,d}.info ; animalb{1,d}.stabilitym];
            idx=(find(conseq==1));
            CP=[animalb{1,d}.info(idx) ; animalb{1,d}.stabilitym(idx)];
            idx=find(conseq==0);
            CN=[animalb{1,d}.info(idx) ; animalb{1,d}.stabilitym(idx)];
            
            
            ControlAll=[ControlAll CA];
            ControlPlace=[ControlPlace CP];
            ControlNonplace=[ControlNonplace CN];
        end
    end
end
%% Philpot
names=Philpotnames;
for n=1:length(names)
    load(char(names(n)));
    for d=1:length(animalb)
        if isempty(animalb{1,d}) == 0
            conseq=animalb{1,d}.conseq;
            PA=[animalb{1,d}.info ; animalb{1,d}.stabilitym];
            idx=(find(conseq==1));
            PP=[animalb{1,d}.info(idx) ; animalb{1,d}.stabilitym(idx)];
            idx=find(conseq==0);
            PN=[animalb{1,d}.info(idx) ; animalb{1,d}.stabilitym(idx)];
            
            
            PhilAll=[PhilAll PA];
            PhilPlace=[PhilPlace PP];
            PhilNonplace=[PhilNonplace PN];
        end
    end
end
%%

%% Plot stability (y) against information content (x)
binnumber=50;
colors = parula(100);
paruladark=colors(1,:);
xlow=0;
xhigh=5;
ylow=-0.5;
yhigh=2;

subplot(3,2,1)
scatterhist(ControlAll(1,:),ControlAll(2,:),'Location','SouthEast','Direction','out','Marker','none')
hold on
h=histogram2(ControlAll(1,:),ControlAll(2,:),binnumber,'DisplayStyle','tile','ShowEmptyBins','on','EdgeAlpha',0);
title('ControlAll')
xlabel('Information content')
ylabel('Stability score')
xlim([xlow xhigh])
ylim([ylow yhigh])
axis square
set(gca,'Color',paruladark)

subplot(3,2,2)
scatterhist(PhilAll(1,:),PhilAll(2,:),'Location','SouthEast','Direction','out','Marker','none')
hold on
h = histogram2(PhilAll(1,:),PhilAll(2,:),binnumber,'DisplayStyle','tile','ShowEmptyBins','on','EdgeAlpha',0);
title('PhilAll')
xlabel('Information content')
ylabel('Stability score')
xlim([xlow xhigh])
ylim([ylow yhigh])
axis square
set(gca,'Color',paruladark)
hold off

subplot(3,2,3)
scatterhist(ControlPlace(1,:),ControlPlace(2,:),'Location','SouthEast','Direction','out','Marker','none')
hold on
h = histogram2(ControlPlace(1,:),ControlPlace(2,:),binnumber,'DisplayStyle','tile','ShowEmptyBins','on','EdgeAlpha',0);
title('ControlPlace')
xlabel('Information content')
ylabel('Stability score')
xlim([xlow xhigh])
ylim([ylow yhigh])
axis square
set(gca,'Color',paruladark)

subplot(3,2,4)
scatterhist(PhilPlace(1,:),PhilPlace(2,:),'Location','SouthEast','Direction','out','Marker','none')
hold on
h = histogram2(PhilPlace(1,:),PhilPlace(2,:),binnumber,'DisplayStyle','tile','ShowEmptyBins','on','EdgeAlpha',0);
title('PhilPlace')
xlabel('Information content')
ylabel('Stability score')
xlim([xlow xhigh])
ylim([ylow yhigh])
axis square
set(gca,'Color',paruladark)

subplot(3,2,5)
scatterhist(ControlNonplace(1,:),ControlNonplace(2,:),'Location','SouthEast','Direction','out','Marker','none')
hold on
h = histogram2(ControlNonplace(1,:),ControlNonplace(2,:),binnumber,'DisplayStyle','tile','ShowEmptyBins','on','EdgeAlpha',0);
title('ControlNonplace')
xlabel('Information content')
ylabel('Stability score')
xlim([xlow xhigh])
ylim([ylow yhigh])
axis square
set(gca,'Color',paruladark)

subplot(3,2,6)
scatterhist(PhilNonplace(1,:),PhilNonplace(2,:),'Location','SouthEast','Direction','out','Marker','none')
hold on
h = histogram2(PhilNonplace(1,:),PhilNonplace(2,:),binnumber,'DisplayStyle','tile','ShowEmptyBins','on','EdgeAlpha',0);
title('PhilNonplace')
xlabel('Information content')
ylabel('Stability score')
xlim([xlow xhigh])
ylim([ylow yhigh])
axis square
set(gca,'Color',paruladark)
