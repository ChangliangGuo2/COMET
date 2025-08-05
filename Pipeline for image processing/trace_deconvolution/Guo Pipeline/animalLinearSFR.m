function [ animal ] = animalLinearSFR( animal, speedThresh, occThresh,binSize )
%ANIMALLINEARSFR Summary of this function goes here
%   Detailed explanation goes here
% sessionNum=1;
% speedThresh=10;
% occThresh=0.05;
% binSize=2;
    for sessionNum=1:length(animal)
        display(['Working on session ' num2str(sessionNum)])
        
        if ~isempty(animal{sessionNum})
            
            animal{sessionNum}.binSize     = binSize;
            animal{sessionNum}.speedThresh = speedThresh;
            animal{sessionNum}.occThresh   = occThresh;
            
            dt = median(diff(animal{sessionNum}.time/1000));            

            position                 = smooth(animal{sessionNum}.pos,ceil(1/dt));%--Smooth the position-----
            dx                       = [0;diff(position)];
            animal{sessionNum}.speed = dx/dt;
            %----For linear track---------------------------------

            animal{sessionNum}.speed = abs(smooth(animal{sessionNum}.speed,ceil(1/dt))).';
            %----Smooth the speed---------------------------------

            %------------Make sure the length of the linear track is 200 cm
            if (~isfield(animal{sessionNum},'trackLength') && max(animal{sessionNum}.pos)>500)
            %--- returns 1 if field is the name of a field of the structure array S. Otherwise, it returns 0.
                animal{sessionNum}.trackLength = 200;
            else %(~isfield(animal{sessionNum},'trackLength'))
                animal{sessionNum}.trackLength = 200;
            end
            %--------------------------------------------------------------
            %% Break position into trials

            % Changed for long track!!!!
            % edges for where to find intersections with the running trace           
            if animal{sessionNum}.trackLength <300
                edge1     = min(animal{sessionNum}.pos)+5;%animal{sessionNum}.trackLength/20;
                edge2     = max(animal{sessionNum}.pos)-15;%animal{sessionNum}.trackLength - animal{sessionNum}.trackLength/20;
            else % 15 is determined by checking the position values. 10 was used. 
                edge1     = 20; %long track: 20
                edge2     = animal{sessionNum}.trackLength - 20;%long track: 20
            end
            %---------------Find the edges of the linear track------------

            tempPos       = position;
            tempPos(isnan(tempPos)) = -100;
            Time_tmp      = animal{sessionNum}.time/1000;
            [t1,p1,idx1]  = polyxpoly(Time_tmp,tempPos,[1 Time_tmp(end)],edge1*[1 1]);
            idx1          = idx1(:,1); %where the position (frame number) intersects with the set edges at the ends
            figure(10000)
                mapshow([1 Time_tmp(end)],edge1*[1 1],'DisplayType','polygon','LineStyle','none')
                hold on;
                mapshow(Time_tmp,tempPos)
                hold on;
                mapshow(t1,p1,'DisplayType','point','Marker','o')               
          

            [t2,p2,idx2]  = polyxpoly(Time_tmp,tempPos,[1 Time_tmp(end)],edge2*[1 1]);
            idx2          = idx2(:,1); %where the position (frame number) intersects with the set edges at the ends
            % combine these two arrays and further process
            figure(10001)
                mapshow([1 Time_tmp(end)],edge2*[1 1],'DisplayType','polygon','LineStyle','none')
                hold on;
                mapshow(Time_tmp,tempPos)
                hold on;
                mapshow(t2,p2,'DisplayType','point','Marker','o')
    
            edgeCrossings      = [];
            edgeCrossings(:,1) = [ones(1,length(idx1)) 2*ones(1,length(idx2))];
            edgeCrossings(:,2) = [idx1' idx2'];
            edgeCrossings      = sortrows(edgeCrossings,2); % sort the matrix by  column 2.
            diffEC             = [diff(edgeCrossings(:,1)); -10]; %add 1 -10 at the end
            %When diffEC is 0, it means the animal is hanging round one
            %edge intead of runing from one edge to the other. When diffEC is
            %1 or -1, it means the animal is running from one edge to the other

            idx1               = find(diffEC==1);
            idx2               = find(diffEC==-1);

            figure(10003)
                yline(edge1)
                hold on;
                mapshow(Time_tmp,tempPos)
                hold on;
                mapshow(t1,p1,'DisplayType','point','Marker','o')
                hold on;
                yline(edge2)
                hold on;               
                mapshow(t2,p2,'DisplayType','point','Marker','square')
                hold on;
                mapshow(Time_tmp(edgeCrossings(idx1,2)),tempPos(edgeCrossings(idx1,2)),'DisplayType','point','Marker','v')
                hold on;
                mapshow(Time_tmp(edgeCrossings(idx2,2)),tempPos(edgeCrossings(idx2,2)),'DisplayType','point','Marker','x')
                hold on;
            % the new idx1 and idx2 are not the ones above. Now they are
            % the index of the matrix (the rows of the matrix below)
            if(idx1(1)>idx2(1))
               idx2 = idx2(2:end); 
            end %continuing to fix trial indexes
            %----Make the left be the start of the trials. So the first
            %trial is when the animal is running from the left(0 cm) to the
            %right (200 cm)

            temp                        = [];
            temp(:,1)                   = edgeCrossings(sort([idx1' idx2']),2);
            temp(:,2)                   = edgeCrossings(sort([idx1' idx2'])+1,2); %temp contains the edges of each trial
            % For example, sort([idx1' idx2']) are 5, 7, 8, 9 and relate to frame
            % number, temp(:,1), are 100, 120, 130, 140, then sort([idx1'
            % idx2'])+1 are 6, 8, 9, 10, so frame number, temp(:,2) are 
            animal{sessionNum}.trialNum = zeros(length(animal{sessionNum}.time),1);
            for i=1:length(temp) 
                animal{sessionNum}.trialNum(temp(i,1):temp(i,2)) = i; %put trial indices on the frames; does not include the portion at the track ends
                if (mod(i,2)==1)
                    plot(Time_tmp(temp(i,1):temp(i,2)),tempPos(temp(i,1):temp(i,2)),'.','Color','blue')
                else
                    plot(Time_tmp(temp(i,1):temp(i,2)),tempPos(temp(i,1):temp(i,2)),'.','Color','green')
                hold on;
                end
            
            end
            %----------the frames from trail 1, 2, 3,......., length(temp)

    %% Find SFR (Spatial firing rate)
            idxSpeed = (animal{sessionNum}.speed>=speedThresh).';

            idx1 = mod(animal{sessionNum}.trialNum,2)==1; %separate intoleft sessions and right sessions
            idx2 = mod(animal{sessionNum}.trialNum+1,2)==1 & animal{sessionNum}.trialNum~=0;

            subs = 1+floor(animal{sessionNum}.pos/binSize);
            subs = subs'; %position vector in units of bins
            animal{sessionNum}.binCenters = (1:ceil(animal{sessionNum}.trackLength/binSize))*binSize-binSize/2;
            
            occ1 = zeros(ceil(animal{sessionNum}.trackLength/binSize),1);
            occ2 = zeros(ceil(animal{sessionNum}.trackLength/binSize),1);           

            %% dealing with broken trials added 2021
        
            ans1=subs(idx1&idxSpeed); %----The bin numbers with frames
            ans1(ans1>100)=100;       %----If the bin number is bigger than 100, then set it to 100.
            temp = accumarray(ans1,dt); %---Caclulate the occupation time in each bin.            
            occ1(1:length(temp)) = temp;            
            occ1 = smoothts(occ1','g',round(30/binSize),5/binSize)'; %gaussian smoothing
            occ1 = occ1';
            animal{sessionNum}.occOdd = occ1;
            occ1(occ1<occThresh) = nan;
            figure(10004)
                subplot(2,1,1);
                heatmap(occ1,'Colormap',parula)


            %same for even trials
            temp = accumarray(subs(idx2&idxSpeed),dt);
            occ2(1:length(temp)) = temp;
            occ2 = smoothts(occ2','g',round(30/binSize),5/binSize)';
            occ2 = occ2';
            animal{sessionNum}.occEven = occ2;
            occ2(occ2<occThresh) = nan;
                subplot(2,1,2);
                h2=heatmap(occ2,'Colormap',parula)


            animal{sessionNum}.FR = nan(length(occ1),size(animal{sessionNum}.firing,2),2);
            temp = zeros(ceil(animal{sessionNum}.trackLength/binSize),size(animal{sessionNum}.firing,2));
            for segNum=1:size(animal{sessionNum}.firing,2)
                
                ans1=subs(idx1&idxSpeed);
                ans1(ans1>100)=100;
                temp2 = accumarray(ans1,animal{sessionNum}.firing(idx1&idxSpeed,segNum)); %accumulate the firing rate of each position bin
                temp(1:length(temp2),segNum) = temp2; % each temp2 might not fill the entire temp map
           
            end
            temp = smoothts(temp','g',round(30/binSize),5/binSize)';
            animal{sessionNum}.FR(:,:,1) = temp./repmat(occ1.',1,size(animal{sessionNum}.firing,2)); %firing rate for each position, divided by occupancy  (otherwise higher rate would be from more time spent in the bin)
            
            %same for even trials:
            temp = zeros(ceil(animal{sessionNum}.trackLength/binSize),size(animal{sessionNum}.firing,2));
            for segNum=1:size(animal{sessionNum}.firing,2)

                temp2 = accumarray(subs(idx2&idxSpeed),animal{sessionNum}.firing(idx2&idxSpeed,segNum));
                temp(1:length(temp2),segNum) = temp2;

            end
            temp = smoothts(temp','g',round(30/binSize),5/binSize)';
            animal{sessionNum}.FR(:,:,2) = temp./repmat(occ2.',1,size(animal{sessionNum}.firing,2)); 
            figure(10005)
                heatmap( animal{sessionNum}.FR(:,1,1)' )
            
            refFR = [animal{sessionNum}.FR(:,:,1); animal{sessionNum}.FR(end:-1:1,:,2)]; %invert the position dimension of half of the trials

            refFR(isnan(refFR)) = 0;
            refFR(refFR==0) = 0.000000000001; %need to get rid of NaNs if you want it to work - redone in the info content calculation
            refFR = refFR./repmat(sum(refFR,1),size(refFR,1),1); %normalize for each cell
            animal{sessionNum}.info = log2(size(refFR,1))+sum(refFR.*log2(refFR),1);
   
            figure(10006)
                plot(real(animal{sessionNum}.info),'-.')
        end
    end
end


