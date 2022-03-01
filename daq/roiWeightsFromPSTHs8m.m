function [holoRequest]=roiWeightsFromPSTHs8m(PSTHs,powers,holoRequest,refPower, samplePeriod)
baselineSubtract=1;
baselinePeriod =1:5;%1:3%
% samplePeriod =8:11;%7:11;%10:14;%
samplePeriod = samplePeriod(1):samplePeriod(2);

sigThreshold = 2;

allPSTHs = [];
allSigCells =[];
figure(1003);clf

powerList = unique(powers);

count=0;
for i=1:numel(powerList)
    p=powerList(i);
    numPassAvg = sum(powers==p);
    x = PSTHs(powers==p,:,:);
    m = squeeze(nanmean(x,1));
    
    x = permute(x,[2 1 3]);
    dataPeriods = x(:,:,samplePeriod);
     sz =size(dataPeriods);
    dataPeriods = reshape(dataPeriods,[sz(1) sz(2)*sz(3)]);
    mData = nanmean(dataPeriods,2);
    
    basePeriods = x(:,:,baselinePeriod);
    sz =size(basePeriods);
    basePeriods = reshape(basePeriods,[sz(1) sz(2)*sz(3)]);
    mBase = nanmean(basePeriods,2);
    sBase = nanstd(basePeriods,[],2)/sqrt(numPassAvg);
    
    threshold = mBase+sigThreshold*sBase;
        threshold2 = mBase-sigThreshold*sBase;

    sigCells = mData>threshold | mData<threshold2;
    allSigCells(:,i)=sigCells;
    
    if baselineSubtract
        b = nanmean(m(:,baselinePeriod),2);
        m = bsxfun(@minus,m,b);
    end
    allPSTHs(:,:,i)=m;
    count=count+1;
    subplot(numel(powerList),2,count)
    plot(m');
    ylabel({['Power : ' num2str(p)]; ['Avg of ' num2str(numPassAvg)]})
        count=count+1;
    subplot(numel(powerList),2,count)
    imagesc(m)
    title([num2str(sum(sigCells)) ' of ' num2str(numel(sigCells)) ' putatively activated'])
end

sigVals = bsxfun(@times,allSigCells,powerList);
sigVals(sigVals==0)=nan;
PowerToUse = min(sigVals,[],2);
if isempty(refPower)
refPower = 0.05;
elseif refPower >0.3
    refPower = refPower/1000;
end

weights = PowerToUse./refPower;

% roiWeights=[];
% for i =1:numel(weights);
%     idx = holoRequest.rois{i};
%     roiWeights(idx) = weights(i);
% end
holoRequest.roiWeights = weights; %in order of targts, not holograms shot.



