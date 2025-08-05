 function [ animal,info ] = animalInfoContent(animal, numShuffles)
    %ANIMALINFOCONTENT Summary of this function goes here
    %   Detailed explanation goes here
 sessionNum       = 1;
 numShuffles      = 500;

    numShuffles   = numShuffles+1;
   for sessionNum = 1:length(animal)
        
       if ~isempty(animal{sessionNum})
            binValueSpread = [];
            
            trialNum = animal{sessionNum}.trialNum;
            
            if mod(max(trialNum),2)==1 %want a even number of trials
                trialNum(trialNum == max(trialNum)) = 0;
            end

            idx1 = mod(trialNum,2)==1; %logical idx for left trials
            idx2 = mod(trialNum+1,2)==1 & trialNum~=0; %right
            
            shift = rand(max(trialNum)/2,numShuffles)*animal{sessionNum}.trackLength*2;
            shift(:,1) = 0; %first info content is real info content
            info  = nan(numShuffles,size(animal{sessionNum}.firing,2));
            for i = 1:numShuffles
                if mod(i,100) == 1
                    display(['Working on info content for session ' num2str(sessionNum) '. ' num2str(i/(numShuffles-1)*100) '% Done'])
                end
                pos = animal{sessionNum}.pos;
                %linearize position for both directions
                pos(idx2) = 2*animal{sessionNum}.trackLength - pos(idx2); %change data to: the mice ran from left to right, and then onwards, then transported back to left beginning
                %------------check how the position is changed!!!!!!!1

                for trial = 1:2:max(trialNum)
                   
                    pos(trialNum == trial | trialNum == (trial+1)) = ... %for each left+right trial
                        mod(pos(trialNum == trial | trialNum == (trial+1)) + ...
                        shift((trial+1)/2,i),2*animal{sessionNum}.trackLength); %using mod, put the shifted position values back to the range of 0~tracklength
     
                end
                
                dt = median(diff(animal{sessionNum}.time/1000));
                
                idxSpeed = (animal{sessionNum}.speed>=animal{sessionNum}.speedThresh);
                
                subs = 1+floor(pos/animal{sessionNum}.binSize);
                subs = subs'; %binned position
                
                occ = zeros(ceil(2*animal{sessionNum}.trackLength/animal{sessionNum}.binSize),1);
                
                temp = accumarray(subs(idxSpeed & trialNum.' ~=0),dt);
                occ(1:length(temp)) = temp;
                occ = smoothts(occ','g',round(30/animal{sessionNum}.binSize),5/animal{sessionNum}.binSize)';
                occ(occ<animal{sessionNum}.occThresh) = nan;
                %calculate SFR as in LinearSFR
                FR = nan(length(occ),size(animal{sessionNum}.firing,2));
                
                temp = zeros(ceil(2*animal{sessionNum}.trackLength/animal{sessionNum}.binSize), ...
                    size(animal{sessionNum}.firing,2));
                for segNum=1:size(animal{sessionNum}.firing,2)
                    temp2 = accumarray(subs(idxSpeed & trialNum.' ~=0),animal{sessionNum}.firing(idxSpeed & trialNum.' ~=0,segNum));
                    temp(1:length(temp2),segNum) = temp2;
                end
                temp = smoothts(temp','g',round(30/animal{sessionNum}.binSize),5/animal{sessionNum}.binSize)';
                FR(:,:) = temp./repmat(occ,1,size(animal{sessionNum}.firing,2));
                
                %spread in bin values
                binValueSpread((end+1):(end+size(FR,1)),:) = FR;
                
                %info content
                refFR = FR;
                refFR(isnan(refFR)) = 0;
                refFR(refFR==0) = 0.000000000001;
                refFR = refFR./repmat(sum(refFR,1),size(refFR,1),1);
                info(i,:) = log2(size(refFR,1))+sum(refFR.*log2(refFR),1);
     
                
            end
            animal{sessionNum}.info = info(1,:); %real info content

            animal{sessionNum}.infoP = sum(repmat(animal{sessionNum}.info,numShuffles-1,1) >= info(2:end,:),1)/(numShuffles-1);
            
            animal{sessionNum}.FR95Per = prctile(binValueSpread,95); % the 95 percentile of the firing rate along the track for each cell
            animal{sessionNum}.FR05Per = prctile(binValueSpread,5);
            animal{sessionNum}.FR50Per = prctile(binValueSpread,50);
            %------------A中已排序的元素被视为第100（0.5/n）、第100（1.5/n），…，第100（[i–0.5]/n）个百分位数。其中i是第i位，表示位置。例如：
            %------------对于诸如{6，3，2，10，1}的五个元素的数据向量，排序的元素{1，2，3，6，10}分别对应于第10、第30、第50、第70和第90百分位数。
            %------------prctile默认使用“exact”线性插值来计算100（0.5/n）和100（[i–0.5]/n）之间的百分比的百分位数
            display(['Percent over 0.95 in session ' num2str(sessionNum) ': ' ...
                num2str(100*sum(animal{sessionNum}.infoP>=0.95)/size(animal{sessionNum}.firing,2)) '%']);
        end
    end
 end