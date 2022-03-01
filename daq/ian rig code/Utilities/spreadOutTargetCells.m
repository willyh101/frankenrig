function cellsOut = spreadOutTargetCells(cellsIn,numCellsPerHolo,holoRequest,verbose)
testvals = unique(cellsIn);
if numel(testvals)==2 && all(testvals==[0 1])
    cellsIn = find(cellsIn);
end
numCells = numel(cellsIn); 
cellIDs = 1:numCells;

Locs = holoRequest.targets;
Locs = Locs(cellsIn,:); 

%% form distance grid
if verbose
disp('Processing Dist Grid')
end

muPerPx = 800/512; %conversion si pixels to microns
muPerOpt = 1/1; %60/14;% %placeholder in case this ratio ever changes

% Locs = hr.targets;
Locs(:,1:2)=Locs(:,1:2).*muPerPx;
Locs(:,3) = Locs(:,3).*muPerOpt;

Dist=[];
for i = 1:size(Locs,1)
    for k=1:size(Locs,1)
        Dist(i,k)=sqrt(sum((Locs(i,:)-Locs(k,:)).^2));
    end
end 

%% make groups
if verbose
disp('First Pass Estimate')
end
numCellsPerHolo;
numGroups = floor(numCells/numCellsPerHolo); 
iterLimit =1e4;

allGroups = nan(numCellsPerHolo,numGroups, iterLimit);
allMinDist = nan([1 iterLimit]);
allMeanDist = nan([1 iterLimit]);

iter=0;
disp('');
while iter<iterLimit
    iter = iter+1;
    
    if mod(iter,1000)==1
        fprintf('.')
    end
    
    randIDs = randperm(numCells);     
    randGroups = reshape(randIDs(1:numCellsPerHolo*numGroups),[numCellsPerHolo numGroups]);
    
    for i=1:numGroups
        tempD = Dist(randGroups(:,i),randGroups(:,i));
        tempD(logical(eye(numCellsPerHolo)))=NaN;
        
        minDist(i) = min(tempD(:));
        meanDist(i) = nanmean(tempD(:));
    end
    
    allGroups(:,:,iter) = randGroups;
    allMinDist(iter) = min(minDist);
    allMeanDist(iter) = mean(meanDist);
    
end

    %% find Best
    
    bestIter = find(allMinDist==max(allMinDist),1);
    
    bestGroup = allGroups(:,:,bestIter);
    bestMin = allMinDist(bestIter);
%     if verbose
        disp(' ')
        disp(['Min Distance First Pass: ' num2str(bestMin)])
%     end
    %% Search a little
    iterLimit2 = 1e4;
    origGroup = bestGroup;
    origMin = bestMin;
    
    iter2 =0;
    while iter2<iterLimit2
        iter2=iter2+1;
        if mod(iter2,1000)==0
            fprintf('.')
        end
        
        newGroup = origGroup;
        
        %transpose
        flip = randperm(numel(newGroup),2);
        hold = newGroup(flip);
        newGroup(flip(2)) = hold(1);
        newGroup(flip(1)) = hold(2);
        
        for i=1:numGroups
            tempD = Dist(newGroup(:,i),newGroup(:,i));
            tempD(logical(eye(numCellsPerHolo)))=NaN;
            
            minDist(i) = min(tempD(:));
            meanDist(i) = nanmean(tempD(:));
        end
        
        thisMin = min(minDist);
        if thisMin>origMin
            if verbose
            disp(['Improvement: min was ' num2str(origMin) ' now ' num2str(thisMin)])
            disp(' ');
            end
            origGroup = newGroup;
            origMin = thisMin;
            iter2=0;
        end
    end
    disp(' ')
    disp(['Final Min Separation: ' num2str(origMin) 'um'])
    
    cellsOut = cellsIn(origGroup(:));
        
    