function [ animal ] = animalStability( animal, numShuffles )
%ANIMALSTABILITY Summary of this function goes here
%   Detailed explanation goes here
% 
%  sessionNum=1;
%  numShuffles=500;

  numShuffles = numShuffles+1;
    for sessionNum=1:length(animal)
        
        if ~isempty(animal{sessionNum})
            binValueSpread = [];
            
            trialNum = animal{sessionNum}.trialNum;
            if mod(max(trialNum),2)==1
                trialNum(trialNum == max(trialNum)) = 0;
            end
            idx1 = mod(trialNum,2)==1;
            idx2 = mod(trialNum+1,2)==1 & trialNum~=0;
            
            idxFirstHalf = trialNum < floor(max(trialNum)/2) & trialNum~=0; %for within session stability
            idxSecondHalf = trialNum > ceil(max(trialNum)/2) & trialNum~=0;
            
            idxOddHalf = (mod(trialNum,4) ==1 | mod(trialNum,4) ==2) & trialNum~=0;
            idxEvenHalf = (mod(trialNum,4) ==3 | mod(trialNum,4) ==0) & trialNum~=0;
            
            shift = rand(max(trialNum)/2,numShuffles)*animal{sessionNum}.trackLength*2;
            shift(:,1) = 0; %first info content is real info content
            stabCorrFirstSecond = nan(numShuffles,size(animal{sessionNum}.firing,2));
            stabCorrOddEven = nan(numShuffles,size(animal{sessionNum}.firing,2));
            for i=1:numShuffles
                if mod(i,100) == 1
                    display(['Working on stability for session ' num2str(sessionNum) '. ' num2str(i/(numShuffles-1)*100) '% Done'])
                end
                pos = animal{sessionNum}.pos;
                %linearize position for both directions
                pos(idx2) = 2*animal{sessionNum}.trackLength - pos(idx2);
                for trial=1:2:max(trialNum) %use OR | to change all the trials. in each shuffle, each trial gets shuffled a different distance
                    pos(trialNum == trial | trialNum == (trial+1)) = ...
                        mod(pos(trialNum == trial | trialNum == (trial+1)) + ...
                        shift((trial+1)/2,i),2*animal{sessionNum}.trackLength);
%                         plot(animal{sessionNum}.pos(trialNum == trial | trialNum == (trial+1)))
%                         hold on
%                         plot(pos(trialNum == trial | trialNum == (trial+1)),'r')
%                         hold off
%                         waitforbuttonpress
                end
                               
                dt = median(diff(animal{sessionNum}.time/1000));
                
                idxSpeed = (animal{sessionNum}.speed>=animal{sessionNum}.speedThresh);            

                subs = 1+floor(pos/animal{sessionNum}.binSize);
                subs = subs';
%% Occupancy:
                occFirstHalf = zeros(ceil(2*animal{sessionNum}.trackLength/animal{sessionNum}.binSize),1);
                
                temp = accumarray(subs(idxSpeed.' & idxFirstHalf),dt); % accumulate how many times the animal were in that position: only for speed above a threshold and left trials
                occFirstHalf(1:length(temp)) = temp;
                occFirstHalf = smoothts(occFirstHalf','g',round(30/animal{sessionNum}.binSize),5/animal{sessionNum}.binSize)';
                occFirstHalf(occFirstHalf<animal{sessionNum}.occThresh) = nan; %exclude places where occupancy was too low
                %similar for second, odd, even
                occSecondHalf = zeros(ceil(2*animal{sessionNum}.trackLength/animal{sessionNum}.binSize),1);
                
                temp = accumarray(subs(idxSpeed.' & idxSecondHalf),dt);
                occSecondHalf(1:length(temp)) = temp;
                occSecondHalf = smoothts(occSecondHalf','g',round(30/animal{sessionNum}.binSize),5/animal{sessionNum}.binSize)';
                occSecondHalf(occSecondHalf<animal{sessionNum}.occThresh) = nan;
                
                occOddHalf = zeros(ceil(2*animal{sessionNum}.trackLength/animal{sessionNum}.binSize),1);
                
                temp = accumarray(subs(idxSpeed.' & idxOddHalf),dt);
                occOddHalf(1:length(temp)) = temp;
                occOddHalf = smoothts(occOddHalf','g',round(30/animal{sessionNum}.binSize),5/animal{sessionNum}.binSize)';
                occOddHalf(occOddHalf<animal{sessionNum}.occThresh) = nan; %?? deleted at first
                
                occEvenHalf = zeros(ceil(2*animal{sessionNum}.trackLength/animal{sessionNum}.binSize),1);
                
                temp = accumarray(subs(idxSpeed.' & idxEvenHalf),dt);
                occEvenHalf(1:length(temp)) = temp;
                occEvenHalf = smoothts(occEvenHalf','g',round(30/animal{sessionNum}.binSize),5/animal{sessionNum}.binSize)';
                occEvenHalf(occEvenHalf<animal{sessionNum}.occThresh) = nan;
                
                
                %calculating firing rates=====================================================
                FRFirstHalf = nan(length(occFirstHalf),size(animal{sessionNum}.firing,2));
                FRSecondHalf = nan(length(occSecondHalf),size(animal{sessionNum}.firing,2));
                
                FROddHalf = nan(length(occOddHalf),size(animal{sessionNum}.firing,2));
                FREvenHalf = nan(length(occEvenHalf),size(animal{sessionNum}.firing,2));
                
                temp = zeros(ceil(2*animal{sessionNum}.trackLength/animal{sessionNum}.binSize), ...
                    size(animal{sessionNum}.firing,2));
                for segNum=1:size(animal{sessionNum}.firing,2)
                    %accumulate firing (for speed > thr, and trial= first half)
                    temp2 = accumarray(subs(idxSpeed.' & idxFirstHalf),animal{sessionNum}.firing(idxSpeed.' & idxFirstHalf,segNum));
                    temp(1:length(temp2),segNum) = temp2;
                end
                temp = smoothts(temp','g',round(30/animal{sessionNum}.binSize),5/animal{sessionNum}.binSize)';
                FRFirstHalf(:,:) = temp./repmat(occFirstHalf,1,size(animal{sessionNum}.firing,2));
                
                temp = zeros(ceil(2*animal{sessionNum}.trackLength/animal{sessionNum}.binSize), ...
                    size(animal{sessionNum}.firing,2));
                for segNum=1:size(animal{sessionNum}.firing,2)
                    temp2 = accumarray(subs(idxSpeed.' & idxSecondHalf),animal{sessionNum}.firing(idxSpeed.' & idxSecondHalf,segNum));
                    temp(1:length(temp2),segNum) = temp2;
                end
                temp = smoothts(temp','g',round(30/animal{sessionNum}.binSize),5/animal{sessionNum}.binSize)';
                FRSecondHalf(:,:) = temp./repmat(occSecondHalf,1,size(animal{sessionNum}.firing,2));
                
                temp = zeros(ceil(2*animal{sessionNum}.trackLength/animal{sessionNum}.binSize), ...
                    size(animal{sessionNum}.firing,2));
                for segNum=1:size(animal{sessionNum}.firing,2)
                    temp2 = accumarray(subs(idxSpeed.' & idxOddHalf),animal{sessionNum}.firing(idxSpeed.' & idxOddHalf,segNum));
                    temp(1:length(temp2),segNum) = temp2;
                end
                temp = smoothts(temp','g',round(30/animal{sessionNum}.binSize),5/animal{sessionNum}.binSize)';
                FROddHalf(:,:) = temp./repmat(occOddHalf,1,size(animal{sessionNum}.firing,2));
                
                temp = zeros(ceil(2*animal{sessionNum}.trackLength/animal{sessionNum}.binSize), ...
                    size(animal{sessionNum}.firing,2));
                for segNum=1:size(animal{sessionNum}.firing,2)
                    temp2 = accumarray(subs(idxSpeed.' & idxEvenHalf),animal{sessionNum}.firing(idxSpeed.' & idxEvenHalf,segNum));
                    temp(1:length(temp2),segNum) = temp2;
                end
                temp = smoothts(temp','g',round(30/animal{sessionNum}.binSize),5/animal{sessionNum}.binSize)';
                FREvenHalf(:,:) = temp./repmat(occEvenHalf,1,size(animal{sessionNum}.firing,2));
                
%                    plot(FRFirstHalf(:,1))
%                    hold on
%                    plot(FRSecondHalf(:,1),'-r')
%                    hold off
%                    waitforbuttonpress
                %info content
                parfor segNum = 1: size(animal{sessionNum}.firing,2)
                    idx = ~isnan(FRFirstHalf(:,segNum)) & ~isnan(FRSecondHalf(:,segNum));
                    idx2 = ~isnan(FROddHalf(:,segNum)) & ~isnan(FREvenHalf(:,segNum));
                    if sum(idx)~=0
%                         stabCorr(i, segNum) = dot(FRFirstHalf(idx,segNum),FRSecondHalf(idx,segNum))/ ...
%                 sqrt(dot(FRFirstHalf(idx,segNum),FRFirstHalf(idx,segNum)) * dot(FRSecondHalf(idx,segNum),FRSecondHalf(idx,segNum)));
                        stabCorrFirstSecond(i, segNum) = corr(FRFirstHalf(idx,segNum),FRSecondHalf(idx,segNum));
                    end
                    if sum(idx2)~=0
                        stabCorrOddEven(i, segNum) = corr(FROddHalf(idx2,segNum),FREvenHalf(idx2,segNum));
                    end
                end
                
                
                
                
            end
            animal{sessionNum}.stabilityCoef = stabCorrFirstSecond(1,:); %real info content
            animal{sessionNum}.stabilityCoefP = sum(repmat(animal{sessionNum}.stabilityCoef,numShuffles-1,1) >= stabCorrFirstSecond(2:end,:),1)/(numShuffles-1);
            
            animal{sessionNum}.stabilityCoef2 = stabCorrOddEven(1,:); %real info content
            animal{sessionNum}.stabilityCoef2P = sum(repmat(animal{sessionNum}.stabilityCoef2,numShuffles-1,1) >= stabCorrOddEven(2:end,:),1)/(numShuffles-1);
            
            
            display(['Percent over 0.95 in session ' num2str(sessionNum) ': ' ...
                num2str(100*sum(animal{sessionNum}.stabilityCoefP>=0.95)/size(animal{sessionNum}.firing,2)) '% | ' ...
                num2str(100*sum(animal{sessionNum}.stabilityCoef2P>=0.95)/size(animal{sessionNum}.firing,2)) '%']);
        end
    end
end


