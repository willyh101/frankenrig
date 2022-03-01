PSTHs = ExpStruct.stimTestDense.PSTHs;
powersToUse = ExpStruct.stimTestDense.powers;

%  PSTHs = PSTHs - mean(PSTHs(:,:,1:5),3);


%%
% highPower = mean(PSTHs(powers==0.1,:,6:end),3);
% figure(4);clf;
% plot(highPower)
% figure(3);clf;
% plot(highPower')
% nanmean(abs(highPower(1,:))./abs(highPower(5,:)))
ErrorCells = mean(mean(isnan(PSTHs(:,:,:)),3),1);
cellsToUse = ErrorCells<0.1;
% cellsToUse = ~isnan(PSTHs(1,:,1));
numCells = sum(cellsToUse);

wholePopResp = nanmean(nanmean(PSTHs(:,:,6:end),3),2);

cellList = find(cellsToUse);
powerList = unique(powersToUse);
pCurve = zeros(numel(cellsToUse),numel(powerList));
pCurveSEM =  zeros(numel(cellsToUse),numel(powerList));

for i = 1:numel(powerList)
    p = powerList(i);
    wholePopMean(i) = nanmean(wholePopResp(powersToUse==p));
end

for i = 1:numel(powerList)
    p = powerList(i);
    theseResp = mean(PSTHs(powersToUse==p,:,6:end),3);
    
    pCurve(:,i) = nanmean(theseResp,1);
    pCurveSEM(:,i) = nanstd(theseResp,1)./sqrt(size(theseResp,1));
end
%%

% [wi hi] = splitPretty(numel(cellsToUse),3,4);
% figure(4);clf
%
% for i=1:numel(cellsToUse);
%     c =i; %cellList(i);
%     subplot(wi,hi,i);
%     e = errorbar(pCurve(c,:),pCurveSEM(c,:));
% end
%

[wi, hi] = splitPretty(numCells,3,4);
figure(4);clf

for i=1:numCells;
    c =cellList(i);
    subplot(wi,hi,i);
    e = errorbar(pCurve(c,:),pCurveSEM(c,:));
    axis tight
end

%%
wholePopResp = nanmean(nanmean(PSTHs(:,:,6:end),3),2);
pidx = zeros(size(powersToUse));
for i = 1:numel(powerList)
    
    p = powerList(i);
    pidx(powersToUse==p) = i;
    wholePopMean(i) = nanmean(wholePopResp(powersToUse==p));
end

figure(6);clf
plot(wholePopResp)
hold on
residualPopResp = wholePopResp' - wholePopMean(pidx);
plot(wholePopResp' - wholePopMean(pidx))

cellMean = squeeze(nanmean(nanmean(PSTHs(:,:,6:end),3),1));

for i = 1:numel(powerList)
    p = powerList(i);
    theseResp = mean(PSTHs(powersToUse==p,:,6:end),3);% - wholePopMean(i);
    
    pCurve(:,i) = nanmean(theseResp,1);
    pCurveSEM(:,i) = nanstd(theseResp,1)./sqrt(size(theseResp,1));
end
figure(4);clf

for i=1:numCells;
    c =cellList(i);
    subplot(wi,hi,i);
    e = errorbar(pCurve(c,:),pCurveSEM(c,:));
    axis tight
    xticklabels([])
    yticklabels([])
end

%% zoom in on particular trace
PSTHsToUse = PSTHs; %ExpStruct.stimTest3.PSTHs; %PSTHs;% - mean(PSTHs(:,:,1:5),3);
PSTHsToUse = powers; %ExpStruct.stimTest3.powers; %powers

cellMean = squeeze(nanmean(nanmean(PSTHsToUse(:,:,6:end),3),1));

wholePopResp = nanmean(nanmean(PSTHsToUse(:,:,6:end),3),2);
pidx = zeros(size(powersToUse));
for i = 1:numel(powerList)
    
    p = powerList(i);
    pidx(powersToUse==p) = i;
    wholePopMean(i) = nanmean(wholePopResp(powersToUse==p));
end

residualPopResp = wholePopResp' - wholePopMean(pidx);


plotIt=0;

if plotIt
figure(5);
end
satVal=[];satValLB=[];Rsquare=[];expoUsed = [];allFits=[];
clear satValExpo satValLBExpo RsquareExpo allFitsExpo
for k = 1:numCells
    if plotIt 
        clf
    end
    cellToPlot = k;
    c =cellList(cellToPlot);
    
    dataToUse = nanmean(PSTHsToUse(:,c,6:end),3) - residualPopResp'*cellMean(c)/(1*nanmean(cellMean));
    
    noLight = nanmean(dataToUse(powersToUse==0));
    dataToUse= dataToUse-noLight;
    
    dataforFit=[];powersforFit=[];
    for i = 1:numel(powerList)
        p = powerList(i);
        theseResp = dataToUse(powersToUse==p);% - wholePopMean(i);
        theseResp(theseResp==max(theseResp)) = [];
        theseResp(theseResp==min(theseResp)) = [];
        
        theseResp(isnan(theseResp))=[];
        
        pCurveExp(i) = nanmean(theseResp,1);
        pCurveSEMExp(i) = nanstd(theseResp,1)./sqrt(size(theseResp,1));
        
        dataforFit = [dataforFit theseResp'];
        powersforFit = [powersforFit ones(size(theseResp'))*p];
    end
    allpCurve(:,k)=pCurveExp;
    allpCurveSEMExp(:,k)=pCurveSEMExp;
    
    if plotIt
        e = errorbar(powerList,pCurveExp,pCurveSEMExp);
        hold on
        scatter(powersToUse,dataToUse,25,residualPopResp,'filled')
        
        %  plot(powerList(2:end),smooth(diff(smooth(pCurveExp,1)),1)*2,'-o')
        refline(0)
    end
    
    expoToTry = [4 6 8];
    
    for expo = 1:numel(expoToTry)
        hillFun = ['(1 / (1 + (-a/x).^' num2str(expoToTry(expo)) ') )*b '];
        [f gof] = fit(powersforFit',dataforFit',hillFun,'Start',[0.1 mean(pCurveExp(end-2:end))] );
        
        
        satValExpo(expo)= f.b;
        ci = confint(f,0.95);
        satValLBExpo(expo) = ci(1,2);
        RsquareExpo(expo) = gof.adjrsquare;
        allFitsExpo{expo} = f;
    end
    expo = find(RsquareExpo == max(RsquareExpo),1);
    satVal(k)= satValExpo(expo);
    ci = confint(f,0.95);
    satValLB(k) =  satValLBExpo(expo);
    Rsquare(k) = RsquareExpo(expo);
    expoUsed(k) = expoToTry(expo);
    
    fToUse = allFitsExpo{expo};
    allFits{k} = fToUse;
    
    %     hillFun = ['(1 / (1 + (-a/x).^' num2str(expoToTry(expo)) ') )*b '];
    %     [f gof] = fit(dataforFit',powersforFit',hillFun,'Start',[0.1 mean(pCurveExp(end-2:end))] );
    
    FatSat(k) = satVal(k)*0.8;
    step = find(fToUse(0.001:0.001:max(powerList)*1.5)>FatSat(k),1);
    if isempty(step)
        PowerAtSat(k)=nan;
    else
        PowerAtSat(k) = step/1000;
    end
    
    
    if plotIt
        plot(allFitsExpo{expo})
        %     pCurveSEMExp(1)
        title(num2str(k))
    end
    %  pause
    fprintf('.')
    if mod(k,50)==0
        disp(num2str(k))
    end
    
end
disp('done')
%% disp
sum(satValLB>0)
figure(6);clf

plot(Rsquare(satValLB>0),satVal(satValLB>0),'o');
cellsThatFit =  satValLB>0 & Rsquare>0.25;
% cellsThatFit =  satValLB>0 & Rsquare>0.25 & Rsquare<0.5;

disp(['Overal R2: ' num2str(mean(Rsquare)) ' R2 used: ' num2str(mean(Rsquare(cellsThatFit)))])
sum(cellsThatFit)


[wi hi] =splitPretty(sum(cellsThatFit),3,4);

figure(4);clf

smallCellList =cellList(cellsThatFit);
cellsThatFitList = find(cellsThatFit);

for i=1:sum(cellsThatFit);
    %     c =smallCellList(i);
    c2 = cellsThatFitList(i);
    subplot(wi,hi,i);
    e = errorbar(powerList,allpCurve(:,c2),allpCurveSEMExp(:,c2));
    axis tight
    hold on
    pwrLine = linspace(min(powerList),max(powerList),100); 
    plot(pwrLine,allFits{c2}(pwrLine),'r');
    plot(PowerAtSat(c2),FatSat(c2),'o')
    hold off
    legend off
    box off
    xlabel ''
    ylabel ''
    drawnow
end

