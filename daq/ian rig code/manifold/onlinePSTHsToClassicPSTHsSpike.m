function [PSTHs, spikes] = onlinePSTHsToClassicPSTHsSpike(onlinePSTHs,onlineEpoch,ExpStruct)
FR = 6.36;

swpStart = ExpStruct.EpochEnterSweep{onlineEpoch};
swpEnd = ExpStruct.EpochEnterSweep{onlineEpoch+1}-1;
swps = swpStart:swpEnd-1;

HRnum = ExpStruct.Holo.Sweeps_holoRequestNumber(swpStart);
HR = ExpStruct.Holo.holoRequests{HRnum};
BLFST = HR.bigListOfFirstStimTimes;

outIDs = ExpStruct.outID(swps);
spikeList = [0 HR.holoStimParams.pulseList];

spikes = outIDs;
uOutIDs = unique(outIDs);
for i = 1:numel(uOutIDs)
    if uOutIDs(i)>0
        spikes(outIDs==uOutIDs(i))=spikeList(i);
    elseif uOutIDs(i)==0
        spikes(outIDs==0) = 0;
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

sweepRunSpeed = cellfun(@(x) mean(computeSpeed(x(:,1))),ExpStruct.digitalSweeps);
RunSpeed = sweepRunSpeed(swps);

if ~isfield(ExpStruct,'runThreshold')
    ExpStruct.runThreshold=inf;
end
runThreshold = ExpStruct.runThreshold;
spikes(RunSpeed>runThreshold)=[];
newPSTHs(RunSpeed>runThreshold,:,:)=[];

PSTHs = newPSTHs;