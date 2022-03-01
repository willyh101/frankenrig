%% get roiWeights from stim test and calculated stimmabiliy

% todo: make this more flexible for loading different holorequests
daqEpochToUse = 4;
startTrial = ExpStruct.EpochEnterSweep{daqEpochToUse}+5; % added this because sometimes the first one is off
hrEpochToUse = ExpStruct.Holo.Sweeps_holoRequestNumber(startTrial);
holoRequest_temp = ExpStruct.Holo.holoRequests{hrEpochToUse};


% out = ExpStruct.sstStimTest;
% out = ExpStruct.pyrStimTest;]
out = ExpStruct.stimTest;
PSTHs=out.PSTHs;
powers = out.powers;

samplePeriod = [8 11];
refPower = 50/1000;
[hr]=roiWeightsFromPSTHs8m(PSTHs,powers,holoRequest_temp,refPower, samplePeriod);
% eligibleCells = find(~isnan(hr.roiWeights) & hr.roiWeights<=2); %find(hr.roiWeights==1);
eligibleCells = find(~isnan(hr.roiWeights) & hr.roiWeights<=2 & hr.roiWeights>=1); %find(hr.roiWeights==1);

ExpStruct.eligibleCells = eligibleCells;
disp([num2str(numel(eligibleCells)) ' stimmable cells']);

%% with power adjustment
holoRequest=hr;
holoRequest.refPower = refPower;
ExpStruct.Holo.holoRequestNumber=ExpStruct.Holo.holoRequestNumber+1;
ExpStruct.Holo.holoRequests{ExpStruct.Holo.holoRequestNumber}=holoRequest;

%% without individual power adj
holoRequest = holoRequest_temp;
ExpStruct.Holo.holoRequestNumber=ExpStruct.Holo.holoRequestNumber+1;
ExpStruct.Holo.holoRequests{ExpStruct.Holo.holoRequestNumber}=holoRequest;

%% make a bunch of single target holos
holosToUse = num2cell(eligibleCells)';
holosToUse
%% make just a few single target holos
holosToUse = [];
nSingleHolos = 30;
cells = randsample(eligibleCells, nSingleHolos);
holosToUse = num2cell(cells)';
holosToUse
%% add a composite holo
% holosToUse{numel(holosToUse)+1} = cells';
holosToUse{numel(holosToUse)+1} = eligibleCells';
holosToUse
%% make a big holo for "1p" stim
holosToUse = [];
nHolos = 5;
nCellsHolo = 50;
for ii=1:nHolos
    cells = randsample(eligibleCells, nCellsHolo);
    holosToUse{numel(holosToUse)+1} = cells';
end
holosToUse