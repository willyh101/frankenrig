out = ExpStruct.stimTest;
PSTHs=out.PSTHs;
powers = out.powers;

% onlinePSTHs = ExpStruct.stimTest.PSTHs;
% [PSTHs, powers] = onlinePSTHsToClassicPSTHs(onlinePSTHs,2,ExpStruct);
% powers = powers(3:46);

% holoRequest.onlinePowerCurve = out;
refPower = 50/1000;
[hr]=roiWeightsFromPSTHs(PSTHs,powers,holoRequest,refPower);
eligibleCells = find(~isnan(hr.roiWeights) & hr.roiWeights<=1); %find(hr.roiWeights==1);
ExpStruct.eligibleCells = eligibleCells;
disp([num2str(numel(eligibleCells)) ' stimmable cells']);

holoRequest=hr;
holoRequest.refPower = refPower;
ExpStruct.Holo.holoRequestNumber=ExpStruct.Holo.holoRequestNumber+1;
ExpStruct.Holo.holoRequests{ExpStruct.Holo.holoRequestNumber}=holoRequest;

visRespCells=[];
%% Vis
PSTHs= ExpStruct.orientationData.PSTHs;
% 
  DAQEpochVis=4;
  [visRespCells vars] = findVisResp(PSTHs,DAQEpochVis,ExpStruct);
  
  %% select cells
  numCellsToSelect = 12;
  
  visAndStim = eligibleCells(ismember(eligibleCells,visRespCells));
  visAndNotStim = eligibleCells(~ismember(eligibleCells,visRespCells));
  
  cellsToUse = [visAndStim' visAndNotStim(randperm(numel(visAndNotStim),numCellsToSelect-numel(visAndStim)))'];
  holosToUse=[];
  for i=1:numCellsToSelect
      holosToUse{i} =  cellsToUse(i);
  end
  
  biggerHolos = spreadOutTargetCells(cellsToUse,10,holoRequest,0); %cellsToCalib; 0; %0; %Put List of Targets to shoot here, or zero to use all available.
  for i=1:numCellsToSelect/10
      holosToUse{end+1} = biggerHolos((i-1)*10+1:i*10);
  end
  
  %% select one at a time v2 (not spread out
  %%form distance grid
muPerPx = 800/512; %conversion si pixels to microns
muPerOpt =65/55; %60/14;% 1/1; %placeholder in case this ratio ever changes

Locs = hr.targets;
Locs(:,1:2)=Locs(:,1:2).*muPerPx;
Locs(:,3) = Locs(:,3).*muPerOpt;

Dist=[];
for i = 1:size(Locs,1)
    for k=1:size(Locs,1)
        Dist(i,k)=sqrt(sum((Locs(i,:)-Locs(k,:)).^2));
    end
end
temp = eye(size(Locs,1));
Dist(logical(temp))=NaN;

%% select Cells
numCellsToSelect = 28;
  

  visAndStim = eligibleCells(ismember(eligibleCells,visRespCells));
  visAndNotStim = eligibleCells(~ismember(eligibleCells,visRespCells));
  
  cellsToUse = [visAndStim' visAndNotStim(randperm(numel(visAndNotStim),numCellsToSelect-numel(visAndStim)))'];
  holosToUse=[];
  
  nHolosEach = 4; %120;
nTargetsPerHolo = 10; %[30]; %[3 5 10 20 33];%[10 3 25];[1 5 10 25];

targetsToUse = cellsToUse; %visCells
nTargets = numel(targetsToUse);

minD = 50;%50; 100;150;
maxD =  inf; 150;300;200;
closeDMin = 50;

c=0;mutalDist = [];
for i=1:numel(nTargetsPerHolo)
    targetsPerHolo = nTargetsPerHolo(i);
    thisSet = nan([nHolosEach targetsPerHolo]);
    thisD = nan([nHolosEach 1]);
    maxReusedCells = floor(targetsPerHolo/2)-1;

    k=1; iter=0;
    while k<=nHolosEach
        nextSet = targetsToUse(randperm(nTargets,targetsPerHolo));
         tempD = Dist(nextSet,nextSet);
         D = nanmean(tempD(:));
         closeD = min(tempD(:));
         
        iter=iter+1;
        if D>maxD || D<minD || closeD < closeDMin
            %don't accept
        elseif ~any(sum(ismember(thisSet,nextSet),2)>=maxReusedCells)...
                && (targetsPerHolo==1 || (D>minD && D<maxD))
            %accept
            thisSet(k,:) = nextSet;
            thisD(k) = D;
            fprintf([' ' num2str(k)])
            if mod(k,25)==0
                disp('.')
            end
            k=k+1;
            
        elseif iter>1e4
            disp('****could not find eligible holo***')
            k=k+1;
            thisSet(k,:) = nan(size(nextSet));
            thisD(k)= nan;
        end
    end
    
    for k=1:nHolosEach
        c=c+1;
        holosToUse{c} = thisSet(k,:);
        mutalDist(c) = thisD(k);
    end
end
  
  
  for i=1:numCellsToSelect
      holosToUse{end+1} =  cellsToUse(i);
  end
  


    
numel(unique([holosToUse{:}]))
 