function [pattern] = complexHoloPattern(target,dura,maxCellPerBin)
%    target = rand([1 300])*5;
%    dura =0.25;
target(target<0)=0;
if size(target,1)~=1
    target = target';
end

freq = 50; %max write rate in Hz
% maxCellPerBin = 20;

nCells = numel(target);
nBins = round(dura*freq);

pattern = zeros([nCells, nBins]);
figure(1);clf
subplot(3,2,6)
histogram(target,100);

ErrList=[];
randOrder = randperm(nCells);
for i=1:nCells; %randOrder;
    val = target(i);
    val = round(val);
    
    est = round(linspace(1,nBins,val));
    pattern(i,est)=1;
    c=0;
    while ~all(sum(pattern)<maxCellPerBin)
        
        c=c+1;
        if c<nBins-1
            pattern(i,est)=0;
            est = round(linspace(1,nBins,val+c));
            est = est+c;%round(linspace(1+c,nBins-c,val));
            est(est>nBins)=[];
            if numel(est)>round(val)
                est(randperm(numel(est),numel(est)-round(val)))=[]; %don't put more spikes than needed
            end
            pattern(i,est)=1;
        else
            fprintf(['Trouble finding space, Target: ' num2str(i) '...']);
            iter=0;
            while iter<1000
                iter=iter+1;
                pattern(i,est)=0;
                est = max(round(rand(1,val).*nBins),1);
                pattern(i,est)=1;
                %                 if numel(unique(est))~=val
                %                     pattern(i,est)=0;
                %                     fprintf([' other problem...']);
                %
                %                     continue
                %                 else
                if numel(unique(est))==val && all(sum(pattern)<maxCellPerBin)
                    fprintf([' Successs!\n']);
                    break
                end
            end
            if ~all(sum(pattern)<maxCellPerBin)
                fprintf('. failed\n')
                %                     disp(['Warning Could not find result! Target: ' num2str(i)])
                pattern(i,est)=0;
                ErrList(end+1)=i;
                break
            end
            
        end
    end
    
    
    subplot(3,2,1);
    imagesc(pattern);
    ylabel('Target')
    xlabel('Bin')
    subplot(3,2,2)
%     hold on
    s = bar(sum(pattern));%,'color',rgb('grey'));
    s.FaceColor = [0.3 0.3 0.3]; %rgb('grey');
    s.BarWidth = 1;
    
    ylabel('Number of Targets Stimmed')
    xlabel('Bin')
%     hold off
    subplot(3,2,3)
%     hold off
    plot(sum(pattern'),target,'o')
    refline(1)
    xlabel('Result N Spikes')
    ylabel('Request N Spikes')
    subplot(3,2,4)
%     hold off
    plot(sum(pattern')-target,'o')
    ylabel('Difference From Delivered to expected')
    xlabel('Target')
    
    subplot(3,2,5)
    minInter=zeros([1 nCells]);meanInter=zeros([1 nCells]);
    for k=1:nCells
        mn = min(diff(find(pattern(k,:))));
        if isempty(mn)
            mn=0;
        end
        minInter(k) = mn;
        mm = mean(diff(find(pattern(k,:))));
        if isempty(mm)
            mm=0;
        end
        meanInter(k)=mm;
        
    end
    fastestFreq = freq./minInter;
    meanFreq = freq./meanInter;
    desiredFreq = target./dura;
    
    p = plot(desiredFreq,fastestFreq,'color','r');
    p.Marker = 'o';
    p.MarkerFaceColor='r';
    p.LineStyle='none';
    refline(1)
    ylabel('Fastest (red)/Mean Freq (black)(Hz)')
    xlabel('Desired Freq (Hz)')
    hold on
    plot(desiredFreq,meanFreq,'o','color','k')
    hold off
    
    
    drawnow
end
disp(['Total Fail Count: ' num2str(numel(ErrList))])