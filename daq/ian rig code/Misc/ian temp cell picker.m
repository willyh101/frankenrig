%% select stim And Vis Cells
nHolosEach = 3;
nTargetsPerHolo = [3 10 20];%[10];% 33];[1 5 10 25];

targetsToUse = stimAndVisCells; stimNotVisCells;
nTargets = numel(targetsToUse);

minD = 100; 100;150;
maxD =  inf;300;200;

c=0;holosToUse=[];mutalDist = [];
for i=1:numel(nTargetsPerHolo)
    targetsPerHolo = nTargetsPerHolo(i);
    thisSet = nan([nHolosEach targetsPerHolo]);
    thisD = nan([nHolosEach 1]);
    maxReusedCells=targetsPerHolo-2;%floor(targetsPerHolo/1.5)-1; %targetsPerHolo-2;%

    k=1; iter=0;
    while k<=nHolosEach
        nextSet = targetsToUse(randperm(nTargets,targetsPerHolo));
         tempD = Dist(nextSet,nextSet);
         D = mean(tempD(:));
        iter=iter+1;
        if ~any(sum(ismember(thisSet,nextSet),2)>=maxReusedCells+1) && (targetsPerHolo==1 || (D>minD && D<maxD))
            thisSet(k,:) = nextSet;
            thisD(k) = D;
            k=k+1;
        elseif iter>1e5
            disp('****could not find eligible holo***')
            k=k+1;
            thisSet(k,:) = nan(size(nextSet));
            thisD(k)= nan;
        end
    end
    
    for k=1:nHolosEach
        c=c+1;
        holosToUse{c} = thisSet(k,:);
        mutalDist(c) = thisD(k);
    end
end
    
numel(unique([holosToUse{:}]))
holosToUse
%%select stim And Vis Cells
nHolosEach = 3;
% nTargetsPerHolo = [3 10 25];%[10];% 33];[1 5 10 25];

targetsToUse = stimNotVisCells;
nTargets = numel(targetsToUse);

minD = 100; 100;150;
maxD =  inf;300;200;

% c=0;holosToUse=[];mutalDist = [];
for i=1:numel(nTargetsPerHolo)
    targetsPerHolo = nTargetsPerHolo(i);
    thisSet = nan([nHolosEach targetsPerHolo]);
    thisD = nan([nHolosEach 1]);
    maxReusedCells=targetsPerHolo-2;%floor(targetsPerHolo/1.5)-1; %targetsPerHolo-2;%

    k=1; iter=0;
    while k<=nHolosEach
        nextSet = targetsToUse(randperm(nTargets,targetsPerHolo));
         tempD = Dist(nextSet,nextSet);
         D = mean(tempD(:));
        iter=iter+1;
        if ~any(sum(ismember(thisSet,nextSet),2)>=maxReusedCells+1) && (targetsPerHolo==1 || (D>minD && D<maxD))
            thisSet(k,:) = nextSet;
            thisD(k) = D;
            k=k+1;
        elseif iter>1e5
            disp('****could not find eligible holo***')
            k=k+1;
            thisSet(k,:) = nan(size(nextSet));
            thisD(k)= nan;
        end
    end
    
    for k=1:nHolosEach
        c=c+1;
        holosToUse{c} = thisSet(k,:);
        mutalDist(c) = thisD(k);
    end
end
    
numel(unique([holosToUse{:}]))
holosToUse
%%select stim And Vis Cells
nHolosEach = 3;
nTargetsPerHolo = [3 10 ];%[10];% 33];[1 5 10 25];

targetsToUse = group1;
nTargets = numel(targetsToUse);

minD = 100; 100;150;
maxD =  inf;300;200;

% c=0;holosToUse=[];mutalDist = [];
for i=1:numel(nTargetsPerHolo)
    targetsPerHolo = nTargetsPerHolo(i);
    thisSet = nan([nHolosEach targetsPerHolo]);
    thisD = nan([nHolosEach 1]);
    maxReusedCells=targetsPerHolo-2;%floor(targetsPerHolo/1.5)-1; %targetsPerHolo-2;%

    k=1; iter=0;
    while k<=nHolosEach
        nextSet = targetsToUse(randperm(nTargets,targetsPerHolo));
         tempD = Dist(nextSet,nextSet);
         D = mean(tempD(:));
        iter=iter+1;
        if ~any(sum(ismember(thisSet,nextSet),2)>=maxReusedCells+1) && (targetsPerHolo==1 || (D>minD && D<maxD))
            thisSet(k,:) = nextSet;
            thisD(k) = D;
            k=k+1;
        elseif iter>1e5
            disp('****could not find eligible holo***')
            k=k+1;
            thisSet(k,:) = nan(size(nextSet));
            thisD(k)= nan;
        end
    end
    
    for k=1:nHolosEach
        c=c+1;
        holosToUse{c} = thisSet(k,:);
        mutalDist(c) = thisD(k);
    end
end
    
numel(unique([holosToUse{:}]))
holosToUse

nHolosEach = 3;
nTargetsPerHolo = [3 10 ];%[10];% 33];[1 5 10 25];

targetsToUse = group2;
nTargets = numel(targetsToUse);

minD = 100; 100;150;
maxD =  inf;300;200;

% c=0;holosToUse=[];mutalDist = [];
for i=1:numel(nTargetsPerHolo)
    targetsPerHolo = nTargetsPerHolo(i);
    thisSet = nan([nHolosEach targetsPerHolo]);
    thisD = nan([nHolosEach 1]);
    maxReusedCells=targetsPerHolo-2;%floor(targetsPerHolo/1.5)-1; %targetsPerHolo-2;%

    k=1; iter=0;
    while k<=nHolosEach
        nextSet = targetsToUse(randperm(nTargets,targetsPerHolo));
         tempD = Dist(nextSet,nextSet);
         D = mean(tempD(:));
        iter=iter+1;
        if ~any(sum(ismember(thisSet,nextSet),2)>=maxReusedCells+1) && (targetsPerHolo==1 || (D>minD && D<maxD))
            thisSet(k,:) = nextSet;
            thisD(k) = D;
            k=k+1;
        elseif iter>1e5
            disp('****could not find eligible holo***')
            k=k+1;
            thisSet(k,:) = nan(size(nextSet));
            thisD(k)= nan;
        end
    end
    
    for k=1:nHolosEach
        c=c+1;
        holosToUse{c} = thisSet(k,:);
        mutalDist(c) = thisD(k);
    end
end
    
numel(unique([holosToUse{:}]))
holosToUse