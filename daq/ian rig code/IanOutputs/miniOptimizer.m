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
    + scoreDistanceRules(x,y)/5 + prioritizeStimmable(x,y)*0.5 ...
    + prioritizeVispResp(x,y)*10 + scoreEnsTuning(x,y) + scoreMeanOSI(x,y);

iterLim = 1e4;

tic
disp('Working on Tuned Holos...')
varsToPass.targetRange = [0 2000];%[0 250];  %distance range based on bounding box (does not work with infs)
varsToPass.ensOSIRange = [0.6 1]; % will optimize to make ensembles in this range

module={};

[finalMatrix, score] = genericGD(startMatrix,iterLim,@updateMatrixFromList,scoreFun,varsToPass,0);
TunedHolos = finalMatrix;

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

tic
disp('Working on unTuned Holos...')
varsToPass.targetRange = [0 2000];%[0 250];  %distance range based on bounding box (does not work with infs)
varsToPass.ensOSIRange = [0 1]; % will optimize to make ensembles in this range

module={};

[finalMatrix, score] = genericGD(startMatrix,iterLim,@updateMatrixFromList,scoreFun,varsToPass,0);
unTunedHolos = finalMatrix;

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
%%
allMatrix = cat(2,TunedHolos,unTunedHolos);
disp(['Using ' num2str(numel(unique(allMatrix(:)))) ' Cells'])
holosToUse = num2cell(allMatrix,1);

holoOptimizer.holosToUse = holosToUse;
holoOptimizer.varsToPass = varsToPass;
holoOptimizer.module =module;

ExpStruct.holoOptimizer = holoOptimizer;