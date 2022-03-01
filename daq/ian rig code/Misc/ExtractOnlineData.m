expEpoch = 10;

swpStrt =ExpStruct.EpochEnterSweep{expEpoch};

try
    swpEnd = ExpStruct.EpochEnterSweep{expEpoch+1}-1;
catch
    swpEnd = numel(ExpStruct.digitalSweeps);
end
swps = swpStrt:swpEnd;

FR = 6.36;


%Run Vector
digitSwp = ExpStruct.digitalSweeps(swps);%% this is where run speed is
digitSwp = cellfun(@(x) x(:,1),digitSwp,'uniformoutput',0); %select line 1 
runVector = cellfun(@(x) computeRunSpeed(x,1/FR),digitSwp,'uniformoutput',0); %convert to cm/s
runVector=cell2mat(runVector');

figure(6);clf
imagesc(runVector)

% idx = floor(1*FR):ceil(2*FR);
% runVal=[];
% for i=1:numel(swps)
%     runVal(i) = mean(runVector(i,idx));
% end

exp.runVal = runVector;

% lowMotionTrials
exp.lowMotionTrials = zeros(size(runVal)); 



% stimID
stimID = ExpStruct.stim_tag(swps);
uniqueStims = unique(stimID);
exp.stimID = stimID; %ExpStruct.outID(swps);

% holo request
HRnum = ExpStruct.Holo.Sweeps_holoRequestNumber(swps(1));
HR = ExpStruct.Holo.holoRequests{HRnum};

% rois
exp.rois = holoRequest.rois;
exp.holoRequest = holoRequest;

try
    bigList =HR.bigListofSequences;
catch
    disp('No bigListofSequences, initializing to no rois')
    bigList=[];
end
uniqueROIs=[];
for i=1:numel(bigList);
    try
    test = bigList{i}; %changed from Exp Struct to holoRequests.biglist b/c a future epoch overwrites this
    catch
        disp('No bigListofSequences, initializing to no rois')
        test=[];
    end
    uniqueROIs{i} = unique(test);

end

stimParam=[];
stimNames =[];
for i =1:numel(uniqueStims)
    ID = uniqueStims(i);
    if ID==0 || ID >=1000
    else
        %         stimNames{i} = ExpStruct.output_names{ID};
%                 stimParam.Seq(i) = sum(diff(ExpStruct.output_patterns{ID}(:,9))>0);
        %         stimParam.numPulse(i) =  sum(diff(ExpStruct.output_patterns{ID}(:,5))>0);
        stimParam.Seq(i) = sum(diff(ExpStruct.stimlog{ID}{1}(:,7))>0);%formerly 9
        stimParam.numPulse(i) =  sum(diff(ExpStruct.stimlog{ID}{1}(:,5))>0);
        if stimParam.Seq(i)~=0
            stimParam.roi{i} = uniqueROIs{stimParam.Seq(i)};
        else
            stimParam.roi{i} = 0;
        end
    end
end
try
    disp('warning holostimparams could be out of order...')
stimParam.Hz = HR.holoStimParams.hzList;
stimParam.numCells = HR.holoStimParams.cellsPerHolo;
stimParam.powers = HR.holoStimParams.powerList;
catch
    disp('no holoStimParams')
end

exp.stimParams = stimParam;

% cell Locations
exp.allCoM = HR.targets(:,1:2);
exp.allDepth = HR.targets(:,3);
%%
invar = msrecv(ExpStruct.SISocket,1)

%% 
out=[];
out.exp = exp;
out.onl = invar;

save(['T:\Ian\OnlineData\OnlineData_' date '.mat'],'out')
disp('saved')
