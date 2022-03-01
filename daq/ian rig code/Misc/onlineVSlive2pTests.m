out = ExpStruct.stimTest; 
PSTHs = out.PSTHs;
powers = out.powers; 

holoRequest.onlinePowerCurve = out;
[hr]=roiWeightsFromPSTHs(PSTHs,powers,holoRequest,50);
eligibleCells = find(~isnan(hr.roiWeights));

out = ExpStruct.stimTest2; 
PSTHs = out.PSTHs;
powers = out.powers; 

holoRequest.onlinePowerCurve = out;
[hr]=roiWeightsFromPSTHs(PSTHs,powers,holoRequest,50);
eligibleCells2 = find(~isnan(hr.roiWeights));

out = ExpStruct.stimTest3; 
PSTHs = out.PSTHs;
powers = out.powers; 

holoRequest.onlinePowerCurve = out;
[hr]=roiWeightsFromPSTHs(PSTHs,powers,holoRequest,50);
eligibleCells3 = find(~isnan(hr.roiWeights));

out = ExpStruct.stimTestRough; 
PSTHs = out.PSTHs;
powers = out.powers; 

holoRequest.onlinePowerCurve = out;
[hr]=roiWeightsFromPSTHs(PSTHs,powers,holoRequest,50);
eligibleCellsR = find(~isnan(hr.roiWeights));

%%
onlineEpoch = 7;
swpStart = ExpStruct.EpochEnterSweep{onlineEpoch};
swpEnd = ExpStruct.EpochEnterSweep{onlineEpoch+1}-1;
swps = swpStart:swpEnd-1;

outIDs = ExpStruct.outID(swps);
powerList = [0:10:100]/1000;

powers = outIDs;
uOutIDs = unique(outIDs);
for i = 1:numel(uOutIDs)
    if uOutIDs(i)>0
        powers(outIDs==uOutIDs(i))=powerList(i);
    elseif uOutIDs(i)==0
        powers(outIDs==0) = 0;
    end
end
        
out = ExpStruct.stimTestOnline; 
PSTHs = out.PSTHs;

HRnum = ExpStruct.Holo.Sweeps_holoRequestNumber(swpStart);
HR = ExpStruct.Holo.holoRequests{6};
BLFST = HR.bigListOfFirstStimTimes;

FR = 6.36;

sz = size(PSTHs);
lastStartFr = round(max(BLFST(:,1))*FR);
minNumFrames = sz(3)-lastStartFr;
firstStartFr = round(min(BLFST(:,1))*FR);

newPSTHs = nan([sz(1) sz(2) minNumFrames+firstStartFr]);
for i =1:sz(2)
    strtFrame =round(BLFST(i,1)*FR);
    if ~isnan(strtFrame)
        newPSTHs(:,i,:) = PSTHs(:,i,strtFrame-firstStartFr+1:strtFrame+minNumFrames);
    end
end

PSTHs = newPSTHs;
% powers = out.powers; 

[hr]=roiWeightsFromPSTHs(PSTHs,powers,holoRequest,50);
eligibleCellsO = find(~isnan(hr.roiWeights));