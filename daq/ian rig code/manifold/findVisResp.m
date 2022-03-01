function [visRespCells vars] = findVisResp(PSTHs,DAQEpoch,ExpStruct)

baselineSubtract = 1;
subtractNull = 0;
runThreshold = 1; %above this excluded
FR = 6.36; %frame Rate in Hz
pAlpha = 0.05;

visLog = cellfun(@(x) length(find(diff(x(:,4))>0)),ExpStruct.digitalSweeps);
startTrial = ExpStruct.EpochEnterSweep{DAQEpoch};

try
    endTrial = ExpStruct.EpochEnterSweep{DAQEpoch+1}-1;
    visLog = visLog(startTrial:endTrial-1);
catch
    endTrial = ExpStruct.sweep_counter;
    % visLog = visLog(startTrial:end);
    visLog = visLog(startTrial:end-1);
end

sweepRunSpeed = cellfun(@(x) mean(computeSpeed(x(:,1))),ExpStruct.digitalSweeps);
visRunSpeed = sweepRunSpeed(startTrial:endTrial-1);

% visLog(visRunSpeed>runThreshold)=[];
% PSTHs(visRunSpeed>runThreshold,:,:)=[];

startTime = round(1 * FR);
endTime = min(round(startTime + 1.5*FR),size(PSTHs,3));

if baselineSubtract
    trialwise_resp = nanmean(PSTHs(:,:,startTime:endTime),3)-nanmean(PSTHs(:,:,1:startTime),3);
else
    trialwise_resp = nanmean(PSTHs(:,:,startTime:endTime),3);
end
if subtractNull
    nullResp = nanmean(nanmean(PSTHs(visLog==min(visLog),:,startTime:endTime),3),1);
    trialwise_resp=bsxfun(@minus,trialwise_resp,nullResp);%trialwise_resp-nullResp;
end
p_vis_resp = [];
for icell = 1:size(trialwise_resp,2)
    this_cell_vals = trialwise_resp(:,icell);
    p_vis_resp(icell) = anova1(this_cell_vals, visLog, 'off');
end

pVals = p_vis_resp;
disp(['There are ' num2str(sum(p_vis_resp<pAlpha)) ' Vis Resp Cells (p<' num2str(pAlpha) ') out of ' num2str(numel(p_vis_resp)) '. aka ' num2str(sum(p_vis_resp<pAlpha)/numel(p_vis_resp)*100,3) '%']);

visRespCells = find(p_vis_resp<pAlpha);
%% Find Responses
vises = unique(visLog);
meanVals=[];
for i=1:numel(vises)
    v =vises(i);
    meanVals(:,i)=mean(trialwise_resp(visLog==v,:));
    eachVal{i} = trialwise_resp(visLog==v,:);
    semVals(:,i)=std(trialwise_resp(visLog==v,:))./sqrt(sum(visLog==v));
end

[maxVisStimVal maxVisStim] = max(meanVals');
figure(10);clf
histogram(maxVisStim(visRespCells),100);
xticks(1:9)
xticklabels({'nan' '0' '45' '90' '135' '180' '225' '270' '315'})

title('All vis Cells')


%% Save some vis variables
vars.visRespCells = visRespCells;
vars.trialwise_resp = trialwise_resp;
vars.visLog = visLog;
vars.meanVals = meanVals;
vars.eachVals = eachVal;
vars.semVals = semVals;

