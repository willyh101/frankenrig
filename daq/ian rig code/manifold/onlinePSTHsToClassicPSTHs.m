function [PSTHs, powers] = onlinePSTHsToClassicPSTHs(onlinePSTHs,onlineEpoch,ExpStruct)
FR = 6.36;

swpStart = ExpStruct.EpochEnterSweep{onlineEpoch};
swpEnd = ExpStruct.EpochEnterSweep{onlineEpoch+1}-1;
swps = swpStart:swpEnd-1;

HRnum = ExpStruct.Holo.Sweeps_holoRequestNumber(swpStart);
try
HR = ExpStruct.Holo.holoRequests{HRnum};
BLFST = HR.bigListOfFirstStimTimes;
catch
    disp('HR error bypassing')
    HR = ExpStruct.Holo.holoRequests{HRnum-1};
BLFST = HR.bigListOfFirstStimTimes;
end
outIDs = ExpStruct.outID(swps);
powerList = [0 HR.holoStimParams.powerList{1}];

powers = outIDs;
uOutIDs = unique(outIDs);
for i = 1:numel(uOutIDs)
    if uOutIDs(i)>0
        powers(outIDs==uOutIDs(i))=powerList(i);
    elseif uOutIDs(i)==0
        powers(outIDs==0) = 0;
    end
end

sz = size(onlinePSTHs);
lastStartFr = round(max(BLFST(:,1))*FR);
minNumFrames = sz(3)-lastStartFr;
firstStartFr = round(min(BLFST(:,1))*FR);

newPSTHs = nan([sz(1) sz(2) minNumFrames+firstStartFr]);
for i =1:sz(2)
    strtFrame =round(BLFST(i,1)*FR);
    if ~isnan(strtFrame)
        newPSTHs(:,i,:) = onlinePSTHs(:,i,strtFrame-firstStartFr+1:strtFrame+minNumFrames);
    end
end

PSTHs = newPSTHs;