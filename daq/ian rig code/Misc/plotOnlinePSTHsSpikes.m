function plotOnlinePSTHsSpikes(PSTHs,spikes,savePath,saveName)
baselineSubtract=1;
baselinePeriod =1:5;
samplePeriod =8:19;

sigThreshold = 2;

if nargin<4
    saveThis=0;
    savePath='';
    saveName='';
    disp('Not Going To Save')
else
    saveThis =1;
end
figure(1004);clf

spikeList = unique(spikes);

count=0;
for i=1:numel(spikeList)
    p=spikeList(i);
    numPassAvg = sum(spikes==p);
    x = PSTHs(spikes==p,:,:);
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
    
    sigCells = mData>threshold;
    
    if baselineSubtract
        b = nanmean(m(:,baselinePeriod),2);
        m = bsxfun(@minus,m,b);
    end
    count=count+1;
    subplot(numel(spikeList),2,count)
    plot(m');

    ylabel({['Num Spikes : ' num2str(p)]; ['Avg of ' num2str(numPassAvg)]})
        count=count+1;
    subplot(numel(spikeList),2,count)
    imagesc(m)
    title([num2str(sum(sigCells)) ' of ' num2str(numel(sigCells)) ' putatively activated'])
    
end

%%%%%%%%%%Rescale axes to match
naxes = length(get(gcf,'Children'));
ymin = 999999; ymax = -ymin;
cmin = 999999; cmax = -cmin;
for k=1:naxes
    if(mod(k,2))
        subplot(naxes/2,2,k)
        ylim = get(gca,'ylim');
        ymin = min([ymin,ylim]);
        ymax = max([ymax,ylim]);
    else
        subplot(naxes/2,2,k)
        clim = get(gca,'clim');
        cmin = min([cmin,clim]);
        cmax = max([cmax,clim]);
    end
end
ylim = [ymin ymax];
clim = [cmin cmax];
for k=1:naxes
    if(mod(k,2))
        subplot(naxes/2,2,k)
        set(gca,'ylim',ylim)
    else
        subplot(naxes/2,2,k)
        set(gca,'clim',clim)
    end
end
%%%%%%%%%%%

out.PSTHs = PSTHs;
out.powers = spikes;

if saveThis
save(fullfile(savePath,saveName),'out')
end