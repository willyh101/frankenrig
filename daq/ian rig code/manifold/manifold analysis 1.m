in = out.spk;

baselineSubtract =0;
subtractNull =1;

FR = info.FR;

stimID = in.stimID;
in.uniqueStims = unique(stimID); 

outID=[];
for i=1:numel(stimID)
    s= stimID(i); 
    s = find(in.uniqueStims==s);
    outID(i) = in.outputsInfo.OutputOrder(s);
end

holosUsed = unique(cat(2,in.holoRequest.bigListofSequences{:}));
roisUsed = unique(cat(2,in.holoRequest.rois{holosUsed})); 
targetsOfInterest = roisUsed; 
% targetsOfInterest = find(~isnan(in.holoRequest.bigListOfFirstStimTimes(:,1)));
% targetsOfInterest = targetsOfInterest( ~isnan(in.targetedCells(targetsOfInterest)));
cellsOfInterest = in.targetedCells(targetsOfInterest);
targetsOfInterest(isnan(cellsOfInterest)) =[];

cellsOfInterest(isnan(cellsOfInterest)) =[];

uniqueOuts = unique(outID);
spikeList =in.holoRequest.holoStimParams.pulseList;

runThreshold = 6;
requireRun =0; %should this be a running mouse or a not running mouse?
if requireRun
    runTrials = mean(in.runVal'>runThreshold)>0.75 ;
else
    runTrials = mean(in.runVal'<runThreshold)>0.75 ;
end

MCThreshold = 3;
Tvector = in.Tarray{1};
Tvector = sqrt(Tvector(:,:,1).^2+Tvector(:,:,1).^2);
newLowMotionTrials = mean(Tvector<MCThreshold)>0.75;


%% Extract Data

dataToUse = in.dfData;
clear mDat
for i=0:numel(spikeList)
    trialsToUse = outID==i & runTrials & newLowMotionTrials ;%& in.lowMotionTrials; 
    data = dataToUse(cellsOfInterest,:,trialsToUse); 
%     data = allDataNoNP(cellsOfInterest,:,trialsToUse); 

if subtractNull
    trialsNull = outID==0 & runTrials & newLowMotionTrials;
    nullData = dataToUse(cellsOfInterest,:,trialsNull); 
end

    sz = size(data);
    
    if i==0
        startTimes = ones(numel(cellsOfInterest),1)*2.5;
    else
        startTimes = in.holoRequest.bigListOfFirstStimTimes(targetsOfInterest,i);
    end
    
    responseData=[];
    for k = 1:numel(cellsOfInterest); 
       strt = startTimes(k)*FR;
       BLPeriod = max(round(strt-0.8*FR),1):round(strt)-1;
       DTPeriod = round(strt+0.3*FR):min(round(strt+1.8*FR),sz(2));
       
       bData =  squeeze(nanmean(data(k,BLPeriod,:),2));
       tData = squeeze(nanmean(data(k,DTPeriod,:),2)); 
       

       if subtractNull
           nData = squeeze(nanmean(nullData(k,DTPeriod,:),2));
           tData = tData-nanmean(nData);
           bData = bData - nanmean(nData); 
       end
       
       if baselineSubtract
           tData = tData-bData;
       end
       
       responseData(k,:) = tData;
    
       
    end
    mDat{i+1}= responseData;
end
%% Fits


counter = 0;

coeff2 = [];
Rsquare2 = [];
fits=[];
plotFast=0;
plotIt=0;

clear minStim maxStim 


if plotIt
    figure(10);clf;
end

mmDatAll=[];
seDatAll=[];

toPlot =1:size(mDat{1},1); 
for i=1:numel(toPlot);%1:size(mDat{1},1);%ExpStruct.eligibleCellsUsedForSpikeTest';%1:size(mDat{1},1)
    
    if mod(i,20)==0
        fprintf([num2str(i) ' ' ])
    end
    
    c=toPlot(i);
    counter=counter+1;
    tempDat = cellfun(@(x) x(c,:),mDat,'uniformoutput',0);
    mmDat = cellfun(@(x) nanmean(x(c,:)),mDat);
    sdDat = cellfun(@(x) nanstd(x(c,:)),mDat);
    nmDat = cellfun(@(x) numel(x(c,:)),mDat);
    seDat = sdDat./sqrt(nmDat);
    
    mmDatAll(counter,:)=mmDat;
    seDatAll(counter,:)=seDat;
    
    spksUsed = [0 spikeList]; %[ 0 1 3 10 30];
    tDat = [];
    gDat = [];
    for i=1:numel(tempDat);
        tDat = [tDat tempDat{i}];
        gDat = [gDat ones([1 numel(tempDat{i})])*spksUsed(i)];
    end
    
    
    % plotSpread(tDat,'categoryIdx',gDat)
    if plotIt
         figure(10);clf;
        subplot(1,2,1)
        plotSpread(tempDat,[],[]);
        e = errorbar(mmDat,seDat);
        e.LineWidth = 2;
        e.Color = rgb('grey');
        title(targetsOfInterest(c));
        ylabel('Fluorescence')
        xlabel('Stim')
        
        
        subplot(1,2,2);
        scatter(gDat,tDat)
    end
    try
        excl = isnan(tDat) | isnan(gDat);
        gDat(excl)=[];
        tDat(excl)=[];
%         ft =fittype({'x.^3','x.^2','x'});
                ft =fittype({'x.^2','x'});

        %         [f2 gof2] = fit(gDat',tDat','poly2');
        %         [f2inv gof2] = fit(tDat',gDat','poly2');
        
                [f2 gof2] = fit(gDat',tDat',ft);
                [f2inv gof2] = fit(tDat',gDat',ft);
%         
%         [f2 gof2] = fit(spksUsed',mmDat',ft);
%         [f2inv gof2] = fit(mmDat',spksUsed',ft);


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
        
    end
    if plotIt
        
        ylabel('Fluorescence')
        xlabel('Pulses added')
        zeroint=fits{counter}(0);
        title({['R2: ' num2str(Rsquare2(counter))]; ['0 intercept: ' num2str(zeroint)]})

    end
    
    minStim(counter) = mmDat(1);
    maxStim(counter) = mmDat(end);

    if ~plotFast &  plotIt
        drawnow
        pause
    end
end
disp('.')


%%
figure(11);clf
subplot(1,2,1)
histogram(Rsquare2,[0:0.02:1])
ylabel('count')
xlabel('R2')
subplot(1,2,2)
% plot(log(abs(coeff2)),log(Rsquare2),'o')
plot(coeff2,Rsquare2,'o')

xlabel('Coeff')
ylabel('RSquared')

%% Co Plot All Fits
figure(12);clf
hold on
xlim([0 20])

toPlot = find(Rsquare2>0.5);
for i=1:numel(toPlot)
    idsx = toPlot(i);
    p = plot(spksUsed,mmDatAll(idsx,:));
    p.Color = [rgb('grey') 0.5];
end
r = refline(0);
r.LineStyle = ':';
r.Color= rgb('grey');

grandMean = mean(mmDatAll(toPlot,:));
grandSEM = std(mmDatAll(toPlot,:))./sqrt(numel(toPlot));
e = errorbar(spksUsed,grandMean,grandSEM);
e.LineWidth =2;
e.Color = 'k';


ylabel('dF/F')
xticks(spksUsed)
xlabel('Holographic Pulse Count')

%% If lost CellIDs
cellIDs = unique(cat(1,holoRequests.rois{:})); 
%% vis section
%calc indiv tuning curves from offline data

%% Write section
%display est spikes

%% Zscore Together
visAllData = vis.allData;
szVis = size(visAllData);
visAllData = reshape(visAllData,[szVis(1) szVis(2)*szVis(3)]);
maniAllData = mani.allData;
szMani = size(maniAllData);
maniAllData = reshape(maniAllData,[szMani(1) szMani(2)*szMani(3)]);


bothAllData  = cat(2,visAllData,maniAllData);

[bothdfData bothZdfData] = computeDFFwithMovingBaseline(bothAllData);

visZdf = bothZdfData(:,1:szVis(2)*szVis(3));
visZdf = reshape(visZdf,szVis);
vis.zdfData2 = visZdf;
maniZdf = bothZdfData(:,szVis(2)*szVis(3)+1:end);
maniZdf = reshape(maniZdf,szMani);
mani.zdfData2 = maniZdf;

%% PcA
% cellsToUse = targettedCells([ 5 49 84 164 177 195 200 205 217 220 232 248 250 257 271 273]);
% targetList =1:numel(targettedCells);

cellsToUse = targettedCells(mani.CellIDs);%1:1000;targettedCells(~ismember(targetList,mani.CellIDs));%targettedCells(visRespCells); targettedCells(mani.CellIDs);%~ROIinArtifact; targettedCells(eligibleCells);
%  cellsToUse = targettedCells(~ismember(targetList,mani.CellIDs));%1:1000;targettedCells(~ismember(targetList,mani.CellIDs));%targettedCells(visRespCells); targettedCells(mani.CellIDs);%~ROIinArtifact; targettedCells(eligibleCells);
% cellsToUse = 1:numCells;%1:1000;targettedCells(~ismember(targetList,mani.CellIDs));%targettedCells(visRespCells); targettedCells(mani.CellIDs);%~ROIinArtifact; targettedCells(eligibleCells);

cellsToUse(isnan(cellsToUse))=[];

visDat = vis.zdfData;
visDat = visDat(cellsToUse,:,:); 

sz= size(visDat);
visDatAll = reshape(visDat,[sz(1) sz(2)*sz(3)]);

[pcaCoef scr] = pca(visDatAll');
pcaCoef1 = pcaCoef(:,1);

% pcaDat = reshape(scr', [sz(1) sz(2) sz(3)]);
% 
pcaDat = pcaCoef' * visDatAll;
pcaDat = reshape(pcaDat,[sz(1) sz(2) sz(3)]);

%%
DAQEpoch = vis.DAQepoch;
visLog = cellfun(@(x) numel(find(diff(x(:,4))==1)),ExpStruct.digitalSweeps,'uniformoutput',1);
startTrial = ExpStruct.EpochEnterSweep{DAQEpoch};

try
    endTrial = ExpStruct.EpochEnterSweep{DAQEpoch+1}-1;
    visLog = visLog(startTrial:endTrial);
catch
    
    % visLog = visLog(startTrial:end);
    visLog = visLog(startTrial:end);
end
% visLog = vis.visID; 

uVis =unique(visLog);
%  uVis=uVis(2:5);
 uVis = [5 8 9];
% clistToUse = colorMapPicker(numel(uVis),'viridis');
clistToUse = colorMapPicker(max(uVis),'viridis');
clistToUse = colorMapPicker(9,'viridis');
% clistToUse = colorMapPicker(5,'viridis');

figure(17);clf
for i =  1:numel(uVis); %[1 4 5]; 1:numel(uVis)
    v = uVis(i);
    pcPlot = mean(pcaDat(:,:,visLog==v | visLog==v+4),3);
    pcPlot = smoothdata(pcPlot,2);
    p = plot3(pcPlot(1,:),pcPlot(2,:),pcPlot(3,:));
%     p = plot(pcPlot(1,:),pcPlot(2,:));
    p.Color = clistToUse{v};
    p.LineWidth=2;
    hold on
end
xlabel('PC1')
ylabel('PC2')
zlabel('PC3')
view(2)
%%
holoDat = mani.zdfData2;
holoDat = holoDat(cellsToUse,:,:);

sz = size(holoDat);
holoDatAll = reshape(holoDat,[sz(1) sz(2)*sz(3)]);

pcaHoloDat = pcaCoef' * holoDatAll;
pcaHoloDat = reshape(pcaHoloDat,[sz(1) sz(2) sz(3)]);


outID = ExpStruct.outID;
swp = ExpStruct.EpochEnterSweep{mani.DAQepoch}: ExpStruct.EpochEnterSweep{mani.DAQepoch+1}-1;
outID = outID(swp);

% stimID = mani.stimID;

colorOrder = 2:9;%mani.outputsInfo.OutputOrder+1;

uHolo = unique(outID);
uHolo = uHolo(2:9);

% uHolo = uHolo([ 2:5] )
for i =  1:numel(uHolo)
    h = uHolo(i);
    pcPlot = mean(pcaHoloDat(:,:,outID==h | outID==h+4),3);
        pcPlot = smoothdata(pcPlot,2);

    p = plot3(pcPlot(1,:),pcPlot(2,:),pcPlot(3,:));
%     p = plot(pcPlot(1,:),pcPlot(2,:));
    p.Color = clistToUse{colorOrder(h)};
    p.LineStyle=':';
    p.LineWidth =2;
    hold on
end
%% Plot mean responses 
figure(10);clf
subplot(1,3,1)

timeEval = round(FR):round((1+1.3)*FR);

% for i=2:5
% dispVisDat(:,i-1) = mean(mean(visDat(:,timeEval,visLog==i ),3),2);
% dispHoloDat(:,i-1) =  mean(mean(holoDat(:,timeEval,outID==5),3),2);
% end
clear dispVisDat dispHoloDat
dispVisDat(:,1) = mean(mean(visDat(:,timeEval,visLog==1),3),2);
dispHoloDat(:,1) =  mean(mean(holoDat(:,timeEval,outID==0),3),2);
for i=2:5;%2:5
dispVisDat(:,i) = mean(mean(visDat(:,timeEval,visLog==i | visLog==i+4),3),2);
dispHoloDat(:,i) =  mean(mean(holoDat(:,timeEval,outID==i-1 | outID==i+3),3),2);
end

imagesc(dispVisDat)
 caxis([-0.5 1.5])
 colorbar
 ylabel('Cells')
xticks(1:5)
xticklabels([nan 0:45:135])
title('Visually Evoked Fluorescence')
 
subplot(1,3,3)
imagesc(dispHoloDat)
 caxis([-0.50 1.5])
 colorbar
 ylabel('Cells')
xticks(1:5)
xticklabels([nan 0:45:135])
title('Holographic Input')
 
 subplot(1,3,2)
 tempEstSpikes(:,1) = estSpikes(:,1);
  tempEstSpikes(:,2:5) = (estSpikes(:,2:5)+estSpikes(:,6:9))./2;;
  imagesc(tempEstSpikes)
  caxis([0 5])
  colorbar
   ylabel('Cells')
xticks(1:5)
xticklabels([nan 0:45:135])
title('Holographic Evoked Fluorescence')

%% if need to reload some stuff
temp = load(physfile,'ExpStruct')
ExpStruct = temp.ExpStruct
mani=exp
eligibleCells = ExpStruct.eligibleCells

mani.CellIDs = unique(cat(1,mani.holoRequest.rois{:}));