%% get roiWeights from Stimtest;
% load('T:\live2p_out\data.mat')
% 
% out.PSTHs = onlinePSTHs;
% out.onlineTraces = onlineTraces;
% out.onlineTraceLengths = onlineTrialLengths;

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
%% vis (this section not used)
% PSTHs= ExpStruct.orientationDataOld.PSTHs;
% PSTHs= ExpStruct.orientationData.PSTHs;
% 
% DAQEpochVis=3;
% [visRespCells vars] = findVisResp(PSTHs,DAQEpochVis,ExpStruct)

%% 

%Define the ensemble Params 
%make a rand group of ensembles

numCellsPerEns = 10;
numEns = 5;

cellsToUse = ExpStruct.eligibleCells;
varsToPass.cellsToUse = cellsToUse;
holoRequest = ExpStruct.Holo.holoRequests{end}; %set this to the correct holoRequest


%generate the random Matrix
startMatrix =zeros([numCellsPerEns,numEns]);
for i=1:numEns
    startMatrix(:,i) = cellsToUse(randperm(numel(cellsToUse),numCellsPerEns));
end

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
 
%% Tuning Rules
visLog = cellfun(@(x) length(find(diff(x(:,4))>0)),ExpStruct.digitalSweeps);

DAQEpochVis = 2;
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

[maxVisStimVal maxVisStim] = max(meanVals');

calc_osi = @(po, oo) ((po - oo)/(po + oo));
tuningCurves = meanVals(:,2:end);
tuningCurves = bsxfun(@minus, tuningCurves,min(tuningCurves,[],2));

[~, prefs] = max(tuningCurves, [], 2);

sz=size(meanVals);
if sz(2)==13;
    orthos1 = prefs-3;
    orthos1(orthos1<1)=orthos1(orthos1<1)+12;
    
    orthos2 = prefs+3;
    orthos2(orthos2>12)=orthos2(orthos2>12)-12;
elseif sz(2)==9
    orthos1 = prefs-2;
    orthos1(orthos1<1)=orthos1(orthos1<1)+8;
    
    orthos2 = prefs+2;
    orthos2(orthos2>8)=orthos2(orthos2>8)-8;
else
    disp('ERROR UNEXPECTED SIZE')
end

orthos = [orthos1 orthos2];
clear osi_vals
for c=1:numel(prefs)
    po = prefs(c);
    oo = orthos(c,:);
    po_vals = tuningCurves(c, po);
    oo_vals = mean(tuningCurves(c,oo));
    osi_vals(c) = calc_osi(po_vals, oo_vals);
end

varsToPass.p_vis_resp =p_vis_resp;
% varsToPass.OriPref = prefs;
varsToPass.osi_vals = osi_vals;
varsToPass.meanOriVals = meanVals; 


%% params to set
varsToPass.percentSimilar = 50; %how similar can two ensembles b?

varsToPass.minSep = 30; %min separation of any two targets
varsToPass.targetRange = [0 300];  %distance range based on bounding box (does not work with infs)

varsToPass.meanOSIrange = [0.5 1]; % will optimize to use individual cells in this range
varsToPass.ensOSIRange = [0.5 1]; % will optimize to make ensembles in this range
varsToPass.pVisLim = 0.05; % will peanalize pVis higher than this value

%%Where we actually build the scoring function
% tic;
scoreFun = @(x,y) checkForTooSimilar(x,y) + avoidHittingCellsToMuch(x)*2 ...
    + scoreDistanceRules(x,y)/5 + prioritizeStimmable(x,y) ...
    + prioritizeVispResp(x,y)*500 + scoreEnsTuning(x,y) + scoreMeanOSI(x,y);

iterLim = 1e4;

%% run the thing GD

% tic
% [finalMatrix, score] = genericGD(startMatrix,iterLim,@updateMatrixFromList,scoreFun,varsToPass,1);
% % [finalMatrix, score] = genericGDparfor(startMatrix,iterLim,@updateMatrixFromList,scoreFun,varsToPass,1);
% % finalMatrix
% 
% disp(['Final Score: ' num2str(score)]);
% disp(['Too Similar Score ' num2str(checkForTooSimilar(finalMatrix,varsToPass))]);
% disp(['Avoid Hitting Cells Too Much ' num2str(avoidHittingCellsToMuch(finalMatrix))]);
% [scr mxDist] =scoreDistanceRules(finalMatrix,varsToPass);
% disp(['Dist Score ' num2str(scr) ' Max Dists: ' num2str(mxDist)]);
% disp(['Prioritize Stimmable score:  ' num2str(prioritizeStimmable(finalMatrix,varsToPass))]);
% [scr PercentAboveThreshold] =prioritizeVispResp(finalMatrix,varsToPass);
% disp(['Prioritize pVis score: ' num2str(scr) '. ' num2str(PercentAboveThreshold*100,2) '% Failed vis Thresh']);
% [scr ensOSIVals ] = scoreEnsTuning(finalMatrix,varsToPass);
% disp(['Ens OSI Score ' num2str(scr) '. Ens OSIs: ' num2str(ensOSIVals)]);
% [scr meanOSIs] = scoreMeanOSI(finalMatrix,varsToPass);
% disp(['Mean OSI Score ' num2str(scr) '. each Ens mean OSIs: ' num2str(meanOSIs)]);
% 
% holosToUse = num2cell(finalMatrix,1); 
% toc



%% Defaults 


% Close Tuned
tic
disp('Working on Close Tuned...')
varsToPass.targetRange = [0 200];%[0 250];  %distance range based on bounding box (does not work with infs)
varsToPass.ensOSIRange = [0.75 1]; % will optimize to make ensembles in this range

module={};


[finalMatrix, score] = genericGD(startMatrix,iterLim,@updateMatrixFromList,scoreFun,varsToPass,0);
closeTuned = finalMatrix;

disp(['Too Similar Score ' num2str(checkForTooSimilar(finalMatrix,varsToPass))]);
disp(['Avoid Hitting Cells Too Much ' num2str(avoidHittingCellsToMuch(finalMatrix))]);
[scr mxDist] =scoreDistanceRules(finalMatrix,varsToPass);
disp(['Dist Score ' num2str(scr) ' Max Dists: ' num2str(mxDist)]);
disp(['Prioritize Stimmable score:  ' num2str(prioritizeStimmable(finalMatrix,varsToPass))]);
[scr PercentAboveThreshold] =prioritizeVispResp(finalMatrix,varsToPass);
disp(['Prioritize pVis score: ' num2str(scr) '. ' num2str(PercentAboveThreshold*100,2) '% Failed vis Thresh']);
[scr ensOSIVals ] = scoreEnsTuning(finalMatrix,varsToPass);
disp(['Ens OSI Score ' num2str(scr) '. Ens OSIs: ' num2str(ensOSIVals)]);
[scr meanOSIs] = scoreMeanOSI(finalMatrix,varsToPass);
disp(['Mean OSI Score ' num2str(scr) '. each Ens mean OSIs: ' num2str(meanOSIs)]);

moduleScores.inputs = varsToPass;
moduleScores.maxDist = mxDist;
moduleScores.ensOSIs = ensOSIVals;
moduleScores.meanOSis = meanOSIs;
moduleScores.score =score;
module{end+1}=moduleScores;
toc
% Close Untuned
disp('Working on Close UNTuned...')

varsToPass.targetRange = [min(mxDist)-50 min(max(mxDist)+50,300)];  %distance range based on bounding box (does not work with infs)
varsToPass.ensOSIRange = [0 0.25]; % will optimize to make ensembles in this range

[finalMatrix, score] = genericGD(startMatrix,iterLim,@updateMatrixFromList,scoreFun,varsToPass,0);
closeUnTuned = finalMatrix;

disp(['Too Similar Score ' num2str(checkForTooSimilar(finalMatrix,varsToPass))]);
disp(['Avoid Hitting Cells Too Much ' num2str(avoidHittingCellsToMuch(finalMatrix))]);
[scr mxDist] =scoreDistanceRules(finalMatrix,varsToPass);
disp(['Dist Score ' num2str(scr) ' Max Dists: ' num2str(mxDist)]);
disp(['Prioritize Stimmable score:  ' num2str(prioritizeStimmable(finalMatrix,varsToPass))]);
[scr PercentAboveThreshold] =prioritizeVispResp(finalMatrix,varsToPass);
disp(['Prioritize pVis score: ' num2str(scr) '. ' num2str(PercentAboveThreshold*100,2) '% Failed vis Thresh']);
[scr ensOSIVals ] = scoreEnsTuning(finalMatrix,varsToPass);
disp(['Ens OSI Score ' num2str(scr) '. Ens OSIs: ' num2str(ensOSIVals)]);
[scr meanOSIs] = scoreMeanOSI(finalMatrix,varsToPass);
disp(['Mean OSI Score ' num2str(scr) '. each Ens mean OSIs: ' num2str(meanOSIs)]);

moduleScores.inputs = varsToPass;
moduleScores.maxDist = mxDist;
moduleScores.ensOSIs = ensOSIVals;
moduleScores.meanOSis = meanOSIs;
moduleScores.score =score;
module{end+1}=moduleScores;
toc

% Far Tuned
disp('Working on Far Tuned...')
varsToPass.targetRange = [500 2000];  %distance range based on bounding box (does not work with infs)
varsToPass.ensOSIRange = [0.75 1]; % will optimize to make ensembles in this range

[finalMatrix, score] = genericGD(startMatrix,iterLim,@updateMatrixFromList,scoreFun,varsToPass,0);
farTuned = finalMatrix;

disp(['Too Similar Score ' num2str(checkForTooSimilar(finalMatrix,varsToPass))]);
disp(['Avoid Hitting Cells Too Much ' num2str(avoidHittingCellsToMuch(finalMatrix))]);
[scr mxDist] =scoreDistanceRules(finalMatrix,varsToPass);
disp(['Dist Score ' num2str(scr) ' Max Dists: ' num2str(mxDist)]);
disp(['Prioritize Stimmable score:  ' num2str(prioritizeStimmable(finalMatrix,varsToPass))]);
[scr PercentAboveThreshold] =prioritizeVispResp(finalMatrix,varsToPass);
disp(['Prioritize pVis score: ' num2str(scr) '. ' num2str(PercentAboveThreshold*100,2) '% Failed vis Thresh']);
[scr ensOSIVals ] = scoreEnsTuning(finalMatrix,varsToPass);
disp(['Ens OSI Score ' num2str(scr) '. Ens OSIs: ' num2str(ensOSIVals)]);
[scr meanOSIs] = scoreMeanOSI(finalMatrix,varsToPass);
disp(['Mean OSI Score ' num2str(scr) '. each Ens mean OSIs: ' num2str(meanOSIs)]);

moduleScores.inputs = varsToPass;
moduleScores.maxDist = mxDist;
moduleScores.ensOSIs = ensOSIVals;
moduleScores.meanOSis = meanOSIs;
moduleScores.score =score;
module{end+1}=moduleScores;

toc

% Far UnTuned
disp('Working on Far UnTuned...')

varsToPass.targetRange = [min(mxDist)-50 max(mxDist)+50];  %distance range based on bounding box (does not work with infs)
varsToPass.ensOSIRange = [0 0.25]; % will optimize to make ensembles in this range

[finalMatrix, score] = genericGD(startMatrix,iterLim,@updateMatrixFromList,scoreFun,varsToPass,0);
FarUnTuned = finalMatrix;

disp(['Too Similar Score ' num2str(checkForTooSimilar(finalMatrix,varsToPass))]);
disp(['Avoid Hitting Cells Too Much ' num2str(avoidHittingCellsToMuch(finalMatrix))]);
[scr mxDist] =scoreDistanceRules(finalMatrix,varsToPass);
disp(['Dist Score ' num2str(scr) ' Max Dists: ' num2str(mxDist)]);
disp(['Prioritize Stimmable score:  ' num2str(prioritizeStimmable(finalMatrix,varsToPass))]);
[scr PercentAboveThreshold] =prioritizeVispResp(finalMatrix,varsToPass);
disp(['Prioritize pVis score: ' num2str(scr) '. ' num2str(PercentAboveThreshold*100,2) '% Failed vis Thresh']);
[scr ensOSIVals ] = scoreEnsTuning(finalMatrix,varsToPass);
disp(['Ens OSI Score ' num2str(scr) '. Ens OSIs: ' num2str(ensOSIVals)]);
[scr meanOSIs] = scoreMeanOSI(finalMatrix,varsToPass);
disp(['Mean OSI Score ' num2str(scr) '. each Ens mean OSIs: ' num2str(meanOSIs)]);
moduleScores.inputs = varsToPass;
moduleScores.maxDist = mxDist;
moduleScores.ensOSIs = ensOSIVals;
moduleScores.meanOSis = meanOSIs;
moduleScores.score =score;
module{end+1}=moduleScores;
toc


%% Final and Export
allMatrix = cat(2,closeTuned,closeUnTuned,farTuned,FarUnTuned);
holosToUse = num2cell(allMatrix,1);

holoOptimizer.holosToUse = holosToUse;
holoOptimizer.varsToPass = varsToPass;
holoOptimizer.module =module;

ExpStruct.holoOptimizer = holoOptimizer;