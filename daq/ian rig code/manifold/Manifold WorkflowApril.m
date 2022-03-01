out = ExpStruct.stimTest; 
PSTHs = out.PSTHs;
powers = out.powers; 

holoRequest.onlinePowerCurve = out;
[hr]=roiWeightsFromPSTHs(PSTHs,powers,holoRequest,50);
eligibleCells = find(~isnan(hr.roiWeights)); %find(hr.roiWeights==1);
noResponseCells = find(isnan(hr.roiWeights));
eligibleCells = find(hr.roiWeights<=2);

ExpStruct.eligibleCells = eligibleCells;
ExpStruct.noResponseCells = noResponseCells;
numel(eligibleCells)
%% prevent roiWeight (like if you're doing hayley pbal later...)
holoRequest.roiWeights = ones(size(holoRequest.roiWeights));

%%
% hr.roiWeights(hr.roiWeights==min(hr.roiWeights)) = 0.5;
holoRequest.roiWeights = hr.roiWeights;

%%
DAQEpoch = 5;
visEpoch = DAQEpoch;
out = ExpStruct.orientationData3;%ExpStruct.ContrastOnlineData%ExpStruct.OrientationData; %; % Select the correct PSTHs
baselineSubtract = 1;
subtractNull = 0;

visLog = cellfun(@(x) length(find(diff(x(:,4))>0)),ExpStruct.digitalSweeps);
startTrial = ExpStruct.EpochEnterSweep{DAQEpoch};

try
    endTrial = ExpStruct.EpochEnterSweep{DAQEpoch+1}-1;
    visLog = visLog(startTrial+1:endTrial);
catch
    
    % visLog = visLog(startTrial:end);
    visLog = visLog(startTrial+1:end );
end


PSTHs = out.PSTHs;
%PSTHs = PSTHs(visLog(visLog>1),:,:);
FR = 6.36; %frame Rate in Hz

startTime = round(1 * FR)+1;

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
pAlpha = 0.05;

disp(['There are ' num2str(sum(p_vis_resp<pAlpha)) ' Vis Resp Cells (p<' num2str(pAlpha) ') out of ' num2str(numel(p_vis_resp)) '. aka ' num2str(sum(p_vis_resp<pAlpha)/numel(p_vis_resp)*100,3) '%']);

%%
visRespCells = p_vis_resp<pAlpha;
 visAndStimCells =eligibleCells(ismember(eligibleCells,find(visRespCells))); 
disp(['There are ' num2str(numel(visAndStimCells)) ' Vis and Stim Resp Cells '])

numToPick =200;

if numel(visAndStimCells)>numToPick
cellsToCalib = visAndStimCells;
else
    notVisCells = eligibleCells(~ismember(eligibleCells,visAndStimCells));
    nToAdd = numToPick-numel(visAndStimCells);
    
    cellsToCalib = [visAndStimCells' notVisCells(randperm(numel(notVisCells),nToAdd))'];
end
disp(['Will calib ' num2str(numel(cellsToCalib)) ' Cells '])

    %%
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
figure;
subplot(1,2,1);
histogram(maxVisStim,100);
title('All Cells')
subplot(1,2,2)
histogram(maxVisStim(visAndStimCells),100);
title('Vis and Stim Cells')

%% Save some vis variables
vars.visRespCells = visRespCells;
vars.visAndStimCells = visAndStimCells;
vars.trialwise_resp = trialwise_resp;
vars.visLog = visLog;
vars.meanVals = meanVals;
vars.eachVals = eachVal;
vars.semVals = semVals;

ExpStruct.visVars = vars;

%% Spike Curve Parts
PSTHs = ExpStruct.spikeTest.PSTHs;
spikes = ExpStruct.spikeTest.spikes;


% plotOnlinePSTHsSpikes(PSTHs,spikes)

%%
% baselineSubtract=0;
baselinePeriod =1:5;
samplePeriod =8:19;

sigThreshold = 2;

spikeList = unique(spikes);

c=0;
for i=1:numel(spikeList)
    p=spikeList(i);
    numPassAvg = sum(spikes==p);
    x = PSTHs(spikes==p,:,:);
    m = squeeze(nanmean(x,1));
    
    x = permute(x,[2 1 3]);
    dataPeriods = nanmean(x(:,:,samplePeriod),3);
    dataPeriodsB = nanmean(x(:,:,baselinePeriod),3);
    
     if subtractNull
        NullData = PSTHs(spikes==0,:,:);
        NullData = permute(NullData,[2 1 3]);
        dataPeriodsNull = nanmean(NullData(:,:,samplePeriod),3);
        dataPeriods = bsxfun(@minus,dataPeriods,nanmean(dataPeriodsNull,2));
        dataPeriodsB = bsxfun(@minus,dataPeriodsB,nanmean(dataPeriodsNull,2));
    end
    if baselineSubtract
        dataPeriods = dataPeriods- dataPeriodsB;
    end
    mDat{i} = dataPeriods;
    [htst , pVals{i}] = ttest(dataPeriodsB', dataPeriods') ;% for stimmability
    
end

%%

counter = 0;

coeff2 = [];
Rsquare2 = [];
fits=[];
plotFast=1;
plotIt=1;0;

if plotIt
    figure(10);
end

clear minStim maxStim minFluor maxFluor


FitCells =cellsToCalib;% look_right;%cellsToCalib;%visAndStimCells; %1:size(mDat{1},1);  ExpStruct.eligibleCells;
ExpStruct.FitCells = FitCells;
for i=1:numel(FitCells);%1:size(mDat{1},1);%ExpStruct.eligibleCellsUsedForSpikeTest';%1:size(mDat{1},1)
    
    if mod(i,20)==0
        fprintf([num2str(i) ' ' ])
    end
    
    c=FitCells(i);
    counter=counter+1;
    if plotIt
    clf
    end
    tempDat = cellfun(@(x) x(c,:),mDat,'uniformoutput',0);
    mmDat = cellfun(@(x) nanmean(x(c,:)),mDat);
    sdDat = cellfun(@(x) nanstd(x(c,:)),mDat);
    nmDat = cellfun(@(x) numel(x(c,:)),mDat);
    seDat = sdDat./sqrt(nmDat);
    
    spksUsed = unique(spikes); %[ 0 1 3 10 30];
    tDat = [];
    gDat = [];
    for i=1:numel(tempDat);
        tDat = [tDat tempDat{i}];
        gDat = [gDat ones([1 numel(tempDat{i})])*spksUsed(i)];
    end
    
    
    % plotSpread(tDat,'categoryIdx',gDat)
    if plotIt
        subplot(1,2,1)
        plotSpread(tempDat,[],[]);
        e = errorbar(mmDat,seDat);
        e.LineWidth = 2;
        e.Color = rgb('grey');
        title(c);
        ylabel('Fluorescence')
        xlabel('Stim')
        
        
        subplot(1,2,2);
        scatter(gDat,tDat)
    end
    try
        % [f1 gof1] = fit(gDat',tDat','poly1');
        excl = isnan(tDat) | isnan(gDat); 
        gDat(excl)=[];
        tDat(excl)=[];
        ft =fittype({'x.^2','x'});
        %         [f2 gof2] = fit(gDat',tDat','poly2');
        %         [f2inv gof2] = fit(tDat',gDat','poly2');
        [f2 gof2] = fit(gDat',tDat',ft);
        [f2inv gof2] = fit(tDat',gDat',ft);
        if plotIt
            
            hold on
            e = errorbar(spksUsed,mmDat,seDat);
            e.LineWidth = 2;
            e.Color = rgb('grey');
            % p1 = plot(f1);
            % p1.Color = 'r';
            % p1.LineWidth=2;
            p2 = plot(f2);
            p2.Color ='k';
            p2.LineWidth=2;
            % legend('Data','Mean',['Linear R2: ' num2str(gof1.rsquare)], ['Quadratic R2: ' num2str(gof2.rsquare) '. coeff: ' num2str(f2.p1)]);
            % legend('Data','Mean',['Linear R2: ' num2str(gof1.rsquare)], ['Quadratic R2: ' num2str(gof2.rsquare) ]);
        end
        
        coeff2(counter) = f2.a; %f2.p1;
        Rsquare2(counter) = gof2.rsquare;
        fits{counter} = f2inv;
        fits2{counter}= f2;
        
    catch
        coeff2(counter) = nan;
        Rsquare2(counter) =0;
        f = fittype('a*x^2 + b*x +c');
        fits{counter} = cfit(f,0,0,0);
%         disp('.')
    end
    if plotIt
        
        ylabel('Fluorescence')
        xlabel('Pulses added')
        
        zeroint=fits{counter}(0);
        title({['R2: ' num2str(Rsquare2(counter))]; ['0 intercept: ' num2str(zeroint)]})
%         subplot(2,2,3)
%         histogram(dat(:,c))
%         xlabel('Fluorescence')
%         ylabel('Frames')
%         title('Measured Fluorescent Values During Visual Stim')
        
        
%         subplot(2,2,4)
%         histogram(fits{counter}(dat(:,c)),0:1:50)
%         xlabel('Est Pulse')
%         ylabel('Count')
%         title('Estimated Number of Pulses Per Frame')
    end
    
    minStim(counter) = mmDat(1);
    maxStim(counter) = mmDat(end);
%     minFluor(counter) = prctile(dat(:,c),10);
%     maxFluor(counter) = prctile(dat(:,c),90);
    if ~plotFast &  plotIt
        drawnow
        pause
    end
end
disp('.')
%%
figure(3);clf
subplot(1,2,1)
histogram(coeff2,100)
ylabel('count')
xlabel('fit coefficient')
subplot(1,2,2)
% plot(log(abs(coeff2)),log(Rsquare2),'o')
plot(coeff2,Rsquare2,'o')

xlabel('Coeff')
ylabel('RSquared')

%% How Many Cells have both a decent R2 and the data is in the range of what was tuned.

% similarRange = minStim < maxFluor & maxStim > minFluor & maxStim > maxFluor & minFluor < minStim;
% sum(similarRange);
% 
% CellsToUse = similarRange & Rsquare2>0.5;
% disp([num2str(sum(similarRange)) ' Similar Range. ' num2str(sum(CellsToUse)) ' also R2>0.5. our of ' num2str(numel(similarRange))]);

%% est spikes for a given mean val
rThreshold =0.5;0.5;

CellsToWrite = find(Rsquare2>rThreshold);
%optional only write vis Resp cells
CellsToWrite = CellsToWrite(visRespCells(FitCells(CellsToWrite)));

disp([num2str(numel(CellsToWrite)) ' have decent fits']);

estSpikes = zeros([numel(CellsToWrite) size(meanVals,2)]);
for  i=1:numel(CellsToWrite)
    cSp = CellsToWrite(i);
    cOri = FitCells(cSp); %visAndStimCells(cSp); 
    estSpikes(i,:) = max(round(fits{cSp}(meanVals(cOri,:))),0);
end
cellIDs = FitCells(CellsToWrite); %visAndStimCells(CellsToWrite);

subtractNullEst=1; %subtract the added spikes observed in the 0 stance

if subtractNullEst
    estSpikes = bsxfun(@minus,estSpikes,estSpikes(:,1));
    estSpikes = max(estSpikes,0);
end
figure(4);
imagesc(estSpikes)
sumEstSpikes = sum(estSpikes)

ExpStruct.rThreshold = rThreshold;
ExpStruct.CellsWritten = cellIDs; 
ExpStruct.estSpikes = estSpikes;

%% same important things
ExpStruct.Mani.CellsToWrite =CellsToWrite;
ExpStruct.Mani.Rsquare2 = Rsquare2;
ExpStruct.Mani.fits = fits; 
ExpStruct.Mani.fits2 = fits2;
ExpStruct.Mani.coeff2 = coeff2;
ExpStruct.Mani.CellsToFit = FitCells; 
ExpStruct.Mani.CellIDs = cellIDs;
ExpStruct.Mani.subtractNullEst = subtractNullEst;

%%
% OriToRep = 4;
% dura=1;
% [pattern] = complexHoloPattern(estSpikes(:,OriToRep),dura);
% %%
% CellsToUse = randi(36,[1 20]);
% fakeData = randi(30,[numel(CellsToUse) 9]);
% [pattern] = complexHoloPattern(fakeData(:,4),1);

%%
% locations = FrankenScopeRigFile();
% load([locations.HoloRequest_DAQ 'holoRequest.mat']);
% disp('loaded holorequest!')

%%

powerList       = {  };
waitList        = [ ]; %time after stim before next stim ms
hzList          = [  ];
pulseList       = [  ]; %previously 10
holosPerCycle   = [  ]; %groups to interleave
cellsPerHolo    = [  ];
divTotalCells   = [  ]; %divide total number of holos
holoSets        = [  ]; %unique groups of cells per
setlinks        = [  ];

rois={};

dura =1;
c=0;
OrisToUse=6:9; %[1:9]; %[6 7 8 9];

maxCellPerBin =5; %default 20

Seq=[];
for o=1:numel(OrisToUse)
    disp(['computing pattern ' num2str(o)])
    Ori = OrisToUse(o);
   
    [pattern] = complexHoloPattern(estSpikes(:,Ori),dura,maxCellPerBin);
    
    
    thisSeq=[];
    for i=1:size(pattern,2)
        cellInHolo = sort(cellIDs(find(pattern(:,i))));
        if ~any(cellfun(@(x) isequal(x,cellInHolo),rois))
            rois{end+1} = cellInHolo;
            thisSeq(i) = numel(rois);
        else
            thisSeq(i) = find(cellfun(@(x) isequal(x,cellInHolo),rois));
        end
    end
    Seq{o} = thisSeq;
    
    
    startTime=1000; %ms; %typically 500
    pulseDuration=5; %ms Stimulation pulse
    TrigDuration = 5; %ms SLM flip command
    stimFreq= 1; % Shouldn't matter but is used
    
    c=c+1;
    powerList(c)       = { 0.1 };
    waitList(c)        = [ 10 ]; %time after stim before next stim ms
    hzList(c)          = [ 50 ];
    pulseList(c)       = [ 1 ]; %previously 10
    holosPerCycle(c)   = [ 1 ]; %groups to interleave
    cellsPerHolo(c)    = [ 0 ];
    divTotalCells(c)   = [ 0 ]; %divide total number of holos
    holoSets(c)        = [ o ]; %unique groups of cells per
    setlinks(c)        = [ 1 ];


% nHolos = numel(rois); 
repsList(c) = numel(thisSeq);
end

%%populate stim params
holoStimParams=[];
holoStimParams.powerList = powerList;
holoStimParams.waitList = waitList;
holoStimParams.hzList = hzList;
holoStimParams.pulseList = pulseList;
holoStimParams.holosPerCycle = holosPerCycle;
holoStimParams.cellsPerHolo = cellsPerHolo;
holoStimParams.divTotalCells = divTotalCells; %not used in subsequent functions
holoStimParams.holoSets = holoSets;
holoStimParams.setlinks = setlinks;
% holoStimParams.nHolos = nHolos;
holoStimParams.repsList = repsList;
% holoStimParams.totalCells = totalCells;
holoStimParams.cellsToUse = cellIDs;

%the other params
holoStimParams.startTime = startTime;
holoStimParams.pulseDuration = pulseDuration;
holoStimParams.TrigDuration = TrigDuration;
holoStimParams.stimFreq = stimFreq;

ExpStruct.Seq =Seq;
holoRequest.rois = rois;

%%
msocketPrep;
holoRequest = transferHR(holoRequest);

ExpStruct.Holo.holoRequestNumber=ExpStruct.Holo.holoRequestNumber+1;
ExpStruct.Holo.holoRequests{ExpStruct.Holo.holoRequestNumber}=holoRequest;

% make the daq sequences
makeHoloTrigSeqs(Seq, holoStimParams, holoRequest);
% save into exp struct
saveExpStructVars(holoRequest, holoStimParams);

%% Check Sweeps for holo trigs
holoTrigTimes = cellfun(@(x) diff(x(:,1)<1000)==1,sweeps,'uniformoutput',0);
holoTrigCounts = cellfun(@(x) sum(diff(x(:,1)<1000)==1),sweeps,'uniformoutput',1);

sweepRunSpeed = cellfun(@(x) mean(computeSpeed(x(:,1))),ExpStruct.digitalSweeps);

%% Check how we did
DAQEpoch = 10;

out = ExpStruct.ManifoldWrite2;

idLog = ExpStruct.outID;
startTrial = ExpStruct.EpochEnterSweep{DAQEpoch};

try
    endTrial = ExpStruct.EpochEnterSweep{DAQEpoch+1}-1;
    idLog = idLog(startTrial:endTrial);
catch
    endTrial = ExpStruct.sweep_counter-1;
    idLog = idLog(startTrial:end);
end

%  idLog=idLog(2:end); 

thisEpochHoloTrigCounts = holoTrigCounts(startTrial:endTrial);

PSTHs = out.PSTHs;
%PSTHs = PSTHs(visLog(visLog>1),:,:);
FR = 6.36; %frame Rate in Hz

startTime = round(1 * FR)+1;

endTime = min(round(startTime + 1.5*FR),size(PSTHs,3));

if baselineSubtract
    maniResp = nanmean(PSTHs(:,:,startTime:endTime),3)-nanmean(PSTHs(:,:,1:startTime),3);
else
    maniResp = nanmean(PSTHs(:,:,startTime:endTime),3);
end

% 
% for icell = 1:size(maniResp,2)
%     maniRespVals = maniResp(:,icell);
% end


ids = unique(idLog);
ids(ids==0)=[];

exclTrials = ~(thisEpochHoloTrigCounts == 50 | thisEpochHoloTrigCounts == 0);
idLog(exclTrials)=-1;

visStartTrial = ExpStruct.EpochEnterSweep{visEpoch};
visEndTrial = ExpStruct.EpochEnterSweep{visEpoch+1}-1;
visEpochRunSpeed = sweepRunSpeed(visStartTrial:visEndTrial);
meanVisRun = mean(visEpochRunSpeed);
stdVisRun = std(visEpochRunSpeed);

thisEpochRunSpeed = sweepRunSpeed(startTrial:endTrial);

abnormalRun = thisEpochRunSpeed>meanVisRun+2*stdVisRun | thisEpochRunSpeed<meanVisRun-2*stdVisRun;
idLog(abnormalRun)= -2; 

idLog(1)=[];

meanManiVals=[];
for i=1:numel(ids)
    v =ids(i);
    meanManiVals(:,i)=mean(maniResp(idLog==v,:));
    semManiVals(:,i)=std(maniResp(idLog==v,:))./sqrt(sum(idLog==v));
end

figure(101);clf
for i=1:numel(cellIDs)
    clf
    writeCell = cellIDs(i);
    fitID = find(FitCells==writeCell);

    subplot(1,3,1)
    tempDat = cellfun(@(x) x(:,writeCell),eachVal,'uniformoutput',0);
    plotSpread(tempDat)
    hold on
    errorbar(meanVals(writeCell,:),semVals(writeCell,:))
    ylabel('Fluorescence (Baselined)')
    xlabel('Ori')
    title(['Cell: ' num2str(writeCell)])
    
    subplot(1,3,2);
    tempDat = cellfun(@(x) x(writeCell,:),mDat,'uniformoutput',0);
    mmDat = cellfun(@(x) nanmean(x(writeCell,:)),mDat);
    sdDat = cellfun(@(x) nanstd(x(writeCell,:)),mDat);
    nmDat = cellfun(@(x) numel(x(writeCell,:)),mDat);
    seDat = sdDat./sqrt(nmDat);
    spksUsed = unique(spikes);
    
    plotSpread(tempDat,'xValues',spksUsed)
     e = errorbar(spksUsed,mmDat,seDat);
        e.LineWidth = 2;
        e.Color = rgb('grey');

        
%         figure;
        plot(fits2{fitID})
            ylabel('Fluorescence')
        xlabel('Stim')
           title(['RSquared = ' num2str(Rsquare2(fitID))]);

        
%     means = [meanVals(writeCell,OrisToUse) meanManiVals(writeCell,oidx)];
%     sems = [semVals(writeCell,ori) semManiVals(writeCell,oidx)];
subplot(1,3,3)   
e =  errorbar(meanVals(writeCell,OrisToUse),semVals(writeCell,OrisToUse));
   e.Color = [0.5 0.5 0.5] ;
   e.LineWidth = 2;
   hold on
   e =  errorbar(meanManiVals(writeCell,1:numel(OrisToUse)),semManiVals(writeCell,1:numel(OrisToUse)));
   e.Color = [1 0 0] ;
   e.LineWidth =2;
   legend('Vis evoked Fluor','Holo evoked Fluor')
   title(['Est Spikes: ' num2str(estSpikes(i,OrisToUse))])
hold off
pause
end


    
