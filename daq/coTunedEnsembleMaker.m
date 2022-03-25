out = ExpStruct.stimTest;
PSTHs=out.PSTHs;
powers = out.powers;


% holoRequest.onlinePowerCurve = out;
refPower = 50/1000;
[hr]=roiWeightsFromPSTHs(PSTHs,powers,holoRequest,refPower);
eligibleCells = find(~isnan(hr.roiWeights) & hr.roiWeights<=2); %find(hr.roiWeights==1);
ExpStruct.eligibleCells = eligibleCells;
disp([num2str(numel(eligibleCells)) ' stimmable cells']);

holoRequest=hr;
holoRequest.refPower = refPower;
ExpStruct.Holo.holoRequestNumber=ExpStruct.Holo.holoRequestNumber+1;
ExpStruct.Holo.holoRequests{ExpStruct.Holo.holoRequestNumber}=holoRequest;

%% Define the ensemble Params 
cellsToUse = ExpStruct.eligibleCells;
varsToPass.cellsToUse = cellsToUse;
holoRequest = ExpStruct.Holo.holoRequests{end}; %set this to the correct holoRequest


%pick your scoring alg/ paramaters that the scorrers will need
pxToMu = 800/512;
zSpacing = 30; %distance between z planes (assumes constant);
targets = holoRequest.targets; 
uniqueZs = unique(targets(:,3));
for i =1:numel(uniqueZs)
    targets(targets(:,3)==uniqueZs(i),3) = (i-1)*zSpacing;
end
targets(:,[1 2]) = targets(:,[1 2]).*pxToMu;

sz = size(targets,1);
Dist = nan([sz(1),sz(1)]);
for i = 1:sz(1)
    for k = 1:sz(1)
        if i==k
            Dist(i,k) = NaN; 
        elseif i<k
            Dist(i,k) = sqrt(sum((targets(i,:)-targets(k,:)).^2));
        end
    end
end
 varsToPass.Dist = Dist;  
 
 % power req
 varsToPass.roiWeights = holoRequest.roiWeights;
 if any(isnan(holoRequest.roiWeights(varsToPass.cellsToUse)))
     disp('ERROR: Some ROIs do not have weights')
 end
 
 
 sz = size(targets,1);
%  center = [256 256 median(targets(:,3))];
center = [256 256];
center = center*pxToMu;
 DistFromCenter = nan(1,[sz(1)]);
for i = 1:sz(1)
    DistFromCenter(i) = sqrt(sum((targets(i,1:2)-center).^2));
end

varsToPass.DistFromCenter = DistFromCenter;
 
%% Tuning Rules
visLog = cellfun(@(x) length(find(diff(x(:,4))>0)),ExpStruct.digitalSweeps);

DAQEpochVis = 3;
startTrial = ExpStruct.EpochEnterSweep{DAQEpochVis};
try
    endTrial = ExpStruct.EpochEnterSweep{DAQEpochVis+1}-1;
    visLog = visLog(startTrial:endTrial-1);
catch
    % visLog = visLog(startTrial:end);
    visLog = visLog(startTrial:end-1);
end

out = ExpStruct.orientationData;
PSTHs = out.PSTHs;
FR = 6.36;
startTime = round(1 * FR)+1;
endTime = min(round(startTime + 1.5*FR),size(PSTHs,3));
trialwise_resp = nanmean(PSTHs(:,:,startTime:endTime),3)-nanmean(PSTHs(:,:,1:startTime),3);
p_vis_resp = [];
for icell = 1:size(trialwise_resp,2)
    this_cell_vals = trialwise_resp(:,icell);
    p_vis_resp(icell) = anova1(this_cell_vals, visLog, 'off');
end
pAlpha = 0.05;

disp(['There are ' num2str(sum(p_vis_resp<pAlpha)) ' Vis Resp Cells (p<' num2str(pAlpha) ') out of ' num2str(numel(p_vis_resp)) '. aka ' num2str(sum(p_vis_resp<pAlpha)/numel(p_vis_resp)*100,3) '%']);
visCells =find(p_vis_resp<pAlpha); 
% sum(ismember(ExpStruct.eligibleCells,visCells))
disp(['There are ' num2str(sum(ismember(ExpStruct.eligibleCells,visCells))) ' vis and stimmable cells'])
% Find Responses
vises = unique(visLog);
meanVals=[];
for i=1:numel(vises)
    v =vises(i);
    meanVals(:,i)=nanmean(trialwise_resp(visLog==v,:));
    semVals(:,i)=std(trialwise_resp(visLog==v,:))./sqrt(sum(visLog==v));
end

calc_osi = @(po, oo) ((po - oo)/(po + oo));
% tuningCurves = meanVals(:,2:end);
tuningCurves = meanVals;
tuningCurves = bsxfun(@minus, tuningCurves,min(tuningCurves,[],2));

%% get prefs, map dirs, go to ori space
[maxVisStimVal, prefs] = max(meanVals');

sz=size(meanVals);
if sz(2)==13
    oris = [nan 0:30:330];
elseif sz(2)==9
    oris = [nan 0:45:315];
else
    disp('ERROR UNEXPECTED SIZE')
end

prefDirs = idx2dirs(prefs);
prefOris = dirs2oris(prefDirs);
% orthoOris = po2oo(prefOris);
orthoOris = pd2oo(prefDirs);

% calculate OSI
clear osi_vals
for c=1:numel(prefs)
%     po = prefOris(c);
    pd = prefDirs(c);
    oo = orthoOris(c,:);
    po_vals = tuningCurves(c, pd==oris);
    oo_vals = mean(tuningCurves(c,ismember(oris, oo)));
    osi_vals(c) = calc_osi(po_vals, oo_vals);
end

figure(51)
clf
histogram(osi_vals(osi_vals > 0))

varsToPass.p_vis_resp =p_vis_resp;
varsToPass.osi_vals = osi_vals;
varsToPass.meanOriVals = meanVals;
varsToPass.oris = oris;
varsToPass.POs = prefs;

%% view some osis

valsToUse = find((osi_vals > 0.5) & (osi_vals < 1));
theseCells = randsample(valsToUse, 16);
disp([num2str(numel(valsToUse)) ' cells'])

figure(66)
clf
for i=1:16
    thisCell = theseCells(i);
    y = meanVals(thisCell,2:end);
    err = semVals(thisCell,2:end);
    subplot(4,4,i)
    errorbar(oris(2:end),y,err)
    title(['Cell #' num2str(thisCell) ' with OSI = ' num2str(osi_vals(thisCell))])
    xlim([0 max(oris)])
end
    


%% params to set
varsToPass.percentSimilar = 50;
varsToPass.minDistFromCenter = 200;

scoreFun = @(x,y) checkForTooSimilar(x,y) + avoidHittingCellsToMuch(x)*2 ...
    + prioritizeStimmable(x,y) + scoreDistanceFromCenter(x,y) + scoreCoTuned(x,y);

iterLim = 1e4;
%%

osi_threshold = 0.4;

% identify the highest co-tuned group
incl = p_vis_resp < pAlpha;
incl = incl & osi_vals > osi_threshold;

[C,~,ic] = unique(prefOris(incl));
[count,idx] = max(accumarray(ic,1));
oriToUse = C(idx);
disp(['Most common ori is ' num2str(oriToUse) ' with ' num2str(count) ' cells.'])

%% make holos
% first make the biggest ensemble possible
cellsToUse = find(prefOris == oriToUse & incl)';
holosToUse{1} = cellsToUse;

%% then use optimizer to get best cases
numEns = 1;

cellRange = 5:count;

clear allScores allEns
c=0;
for j=cellRange
    c=c+1;
    startMatrix =zeros([j,numEns]);
    for i=1:numEns
        startMatrix(:,i) = cellsToUse(randperm(numel(cellsToUse),j));
    end
    [finalMatrix, score] = genericGD(startMatrix,iterLim,@updateMatrixFromList,scoreFun,varsToPass,0);
    allEns{c} = finalMatrix;
    allScores(c) = score;
    disp([num2str(j) ' cells in ensemble, score: ' num2str(score)])
end

normScores = allScores./sqrt(cellRange);
[~, sidx] = sort(normScores);
% [~, sidx] = sort(allScores);

disp('These are your best ensemble sizes...')
disp(cellRange(sidx))

%%
numEns = 1;
szEns = cellRange(sidx(1:3));
szEns = [szEns 10 3];

clear allScores allEns
c=0;
for j=szEns
    c=c+1;
    startMatrix =zeros([j,numEns]);
    for i=1:numEns
        startMatrix(:,i) = cellsToUse(randperm(numel(cellsToUse),j));
    end
    [finalMatrix, score] = genericGD(startMatrix,iterLim,@updateMatrixFromList,scoreFun,varsToPass,0);
    allEns{c} = finalMatrix;
    allScores(c) = score;
    disp(score)
end

%%
k=numel(holosToUse);
nHolosAdd = numEns*numel(szEns);
for i=1:numel(szEns)
    temp = allEns{i};
    for j=1:numEns
        k=k+1;
        holosToUse{k} = temp(:,j);
    end
end

%%
figure(11)
clf
hold on
scatter(center(2), center(1))
for i=1:numel(holosToUse)
    this_holo = holosToUse{i};
    locs = targets(this_holo,:);
    scatter(locs(:,2), locs(:,1), 'filled')
    xlim([0,512*pxToMu])
    ylim([0,512*pxToMu])
end

%%
ExpStruct.isoOriShown = 90;
ExpStruct.sizes = [10 20 50];
ExpStruct.contrast = 100;
ExpStruct.oris = [nan 90 0];
ExpStruct.gratingLoc = 'center';