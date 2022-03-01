%% Online Data
mDatOnline = mDat;
%% Offline data
mDatOffline = mDat;


%%

figure(10); clf

toPlot =targetsOfInterest; %1:size(mDat{1},1);  ExpStruct.eligibleCells;
for ii=1:numel(toPlot);%1:size(mDat{1},1);%ExpStruct.eligibleCellsUsedForSpikeTest';%1:size(mDat{1},1)
    c=toPlot(ii);
    counter=counter+1;
    clf
    tempDat = cellfun(@(x) x(c,:),mDatOnline,'uniformoutput',0);
    mmDat = cellfun(@(x) nanmean(x(c,:)),mDatOnline);
    sdDat = cellfun(@(x) nanstd(x(c,:)),mDatOnline);
    nmDat = cellfun(@(x) numel(x(c,:)),mDatOnline);
    seDat = sdDat./sqrt(nmDat);
    
    spksUsed = unique(spikes); %[ 0 1 3 10 30];
    tDat = [];
    gDat = [];
    for i=1:numel(tempDat);
        tDat = [tDat tempDat{i}];
        gDat = [gDat ones([1 numel(tempDat{i})])*spksUsed(i)];
    end
    
    
    if plotIt
        subplot(2,2,1)
        plotSpread(tempDat,[],[]);
        e = errorbar(mmDat,seDat);
        e.LineWidth = 2;
        e.Color = rgb('grey');
        title(c);
        ylabel('Fluorescence')
        xlabel('Stim')
        
        
        subplot(2,2,2);
        scatter(gDat,tDat)
    end
    try
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
            p2 = plot(f2);
            p2.Color ='k';
            p2.LineWidth=2;
        end
        
        coeff2(counter) = f2.a;
        Rsquare2(counter) = gof2.rsquare;
        fits{counter} = f2inv;
        
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
    
    c=ii;
    counter=counter+1;
    tempDat = cellfun(@(x) x(c,:),mDatOffline,'uniformoutput',0);
    mmDat = cellfun(@(x) nanmean(x(c,:)),mDatOffline);
    sdDat = cellfun(@(x) nanstd(x(c,:)),mDatOffline);
    nmDat = cellfun(@(x) numel(x(c,:)),mDatOffline);
    seDat = sdDat./sqrt(nmDat);
    
    spksUsed = unique(spikes); %[ 0 1 3 10 30];
    tDat = [];
    gDat = [];
    for i=1:numel(tempDat);
        tDat = [tDat tempDat{i}];
        gDat = [gDat ones([1 numel(tempDat{i})])*spksUsed(i)];
    end
    
    
    if plotIt
        subplot(2,2,3)
        plotSpread(tempDat,[],[]);
        e = errorbar(mmDat,seDat);
        e.LineWidth = 2;
        e.Color = rgb('grey');
        title(targetsOfInterest(c));
        ylabel('Fluorescence')
        xlabel('Stim')
        
        
        subplot(2,2,4);
        scatter(gDat,tDat)
    end
    try
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
            p2 = plot(f2);
            p2.Color ='k';
            p2.LineWidth=2;
        end
        
        coeff2(counter) = f2.a;
        Rsquare2(counter) = gof2.rsquare;
        fits{counter} = f2inv;
        
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
    drawnow
    pause
end