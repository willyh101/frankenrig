%% if desired, load hr
locations = FrankenScopeRigFile();
load([locations.HoloRequest_DAQ 'holoRequest.mat']);

%% Manually set holosToUse
holosToUse{1} = [1:5];
holosToUse{2} = [1 2 3]; 

holosToUse{1} = [2 5 4 9 8  3 13 6 7 1   20 23]; %all dec
holosToUse{2} = [2 5 4 9 8] ; %best 1 plane
holosToUse{1} = [2 5 4 9 8  3 13 6 7 1]; %dec 1 plane

% holosToUse{3} = 3; 
% holosToUse{4} = 4; 
% holosToUse{5} = 5; 
holosToUse{1} = [1:36];
% holosToUse{4} = 33; 
holosToUse{1} = [1:638];
holosToUse{1} = [1:size(holoRequest.targets,1)];

holosToUse{5} = [1:11 32 17 19 60 70 82 84 20 12];

holosToUse{1} = [6 31 48    2 1 4];
holosToUse{1} = [65 66 68 71 86 75   62 63 64];
holosToUse{1} = [131 157 147 151    123 124 127];

holosToUse=[];
holosToUse{1} = [1 2 3];[6 31 48];
holosToUse{2} = [65 66 68];
holosToUse{3} = [2 1 4];
%% 
targets = holoRequest.targets;

limY = [70 320]; %[1 350];
limX = [135 340]; %[40 450];
limZ = [0 0];

targToUse = find( targets(:,1)>=limY(1) & targets(:,1)<=limY(2) ...
    & targets(:,2)>=limX(1) & targets(:,2)<=limX(2) ...
    & targets(:,3)>=limZ(1) & targets(:,3)<=limZ(2) );

holosToUse=[];
holosToUse{1} = targToUse';
holosToUse{2} = targToUse(randperm(numel(targToUse),10))';
% holosToUse{3} = targToUse(randperm(numel(targToUse),20))';
% holosToUse{4} = targToUse(randperm(numel(targToUse),50))';

limZ = [15 15];
targToUse = find( targets(:,1)>=limY(1) & targets(:,1)<=limY(2) ...
    & targets(:,2)>=limX(1) & targets(:,2)<=limX(2) ...
    & targets(:,3)>=limZ(1) & targets(:,3)<=limZ(2) );
holosToUse{3} = targToUse(randperm(numel(targToUse),10))';

limZ = [30 30];
targToUse = find( targets(:,1)>=limY(1) & targets(:,1)<=limY(2) ...
    & targets(:,2)>=limX(1) & targets(:,2)<=limX(2) ...
    & targets(:,3)>=limZ(1) & targets(:,3)<=limZ(2) );
holosToUse{4} = targToUse(randperm(numel(targToUse),10))';

holosToUse{5} = [holosToUse{2:4}];

% holosToUse{1} = targToUse(randperm(numel(targToUse),20))';

%% set params
startTime = 500; %ms;
visStartTime = 1000; %ms %added 4/5/21 -IAO
pulseDuration = 3000; %ms Stimulation pulse 
TrigDuration = 5; %ms SLM flip command
stimFreq= 1; % Shouldn't matter but is used

P1=[0.05 0.1 0.2 0.3]; %[0.001 .005 0.01 .025 .05 .1]; %[.05 .1 .15]; %[0.01 0.025 0.05 .075 0.1 .15]; %List of Powers To Use

numHolos=numel(holosToUse);
% powerList       = { P1 };
% waitList        = [ 6]; %time after stim before next stim ms
% hzList          = [ 30 ];
% pulseList       = [ 10 ]; %previously 10
% holosPerCycle   = [ 1 ]; %groups to interleave
% cellsPerHolo    = [ 3 ];
% divTotalCells   = [ 1 ]; %divide total number of holos 
% holoSets        = [ 1 ]; %unique groups of cells per 
% setlinks        = [ 1 ];
clear powerList waitList hzList pulseList holosPerCycle cellsPerHolo divTotalCells holoSets setlinks
for i=1:numHolos
    powerList{i}       = [P1];
    waitList(i)        = 10; %time after stim before next stim ms
    hzList(i)          = 30;  %30Hz
    pulseList(i)       = 1; %5 is default
    holosPerCycle(i)   = 1; %groups to interleave
    cellsPerHolo(i)    = numel(holosToUse{i});
    divTotalCells(i)   = numHolos; %divide total number of holos
    holoSets(i)        = i; %unique groups of cells per
    setlinks(i)        = 1;
end
cellsToUse =holosToUse;%repmat(holosToUse,[1 numel(stimTypeNumPulse)]);%[holosToUse{:}];% holosToUse;%GroupList;%cellsToShoot;eligibleCells;% 0; %Put List of Targets to shoot here, or zero to use all available.

if ~iscell(cellsToUse) && cellsToUse==0
    totalCells =size(holoRequest.targets,1);
    disp(['Total Cells Detected ' num2str(totalCells)]);
    cellsToUse = 1:totalCells;
elseif iscell(cellsToUse)
    totalCells = numel(unique([cellsToUse{:}]));
    disp(['Using ' num2str(totalCells) ' Cells']);
else
    totalCells = numel(cellsToUse);
    disp(['Using ' num2str(totalCells) ' Cells']);
end

nHolos          = floor(totalCells./divTotalCells./cellsPerHolo); %only make complete holograms
repsList        = divTotalCells./divTotalCells;%floor(nHolos./holosPerCycle); 

SequenceDurations =repsList./hzList.*pulseList+startTime/1000; %in s
maxSeqDur =max(SequenceDurations);
disp(['The Longest Sequence is ' num2str(maxSeqDur,3) 's']);

%%populate stim params
holoStimParams.powerList = powerList;
holoStimParams.waitList = waitList;
holoStimParams.hzList = hzList;
holoStimParams.pulseList = pulseList;
holoStimParams.holosPerCycle = holosPerCycle;
holoStimParams.cellsPerHolo = cellsPerHolo;
holoStimParams.divTotalCells = divTotalCells;
holoStimParams.holoSets = holoSets;
holoStimParams.setlinks = setlinks;
% holoStimParams.nHolos = nHolos;
holoStimParams.repsList = repsList;
holoStimParams.totalCells = totalCells;
holoStimParams.cellsToUse = cellsToUse;

%the other params
holoStimParams.startTime = startTime;
holoStimParams.pulseDuration = pulseDuration;
holoStimParams.TrigDuration = TrigDuration;
holoStimParams.stimFreq = stimFreq;
holoStimParams.visStartTime = visStartTime; %added 4/5/21 -IAO

holoRequest.holoStimParams = holoStimParams;

setKey = {};
rois = {};
for i = 1:numHolos
    rois{i} = cellsToUse{i}; %cellsToUse(i,:);
    setKey{i} = i;
end

holoRequest.rois = rois;

%% run the universal stuff
% talk to holo comp, fill in roi weights if its not there
msocketPrep;
holoRequest = transferHR(holoRequest);

ExpStruct.Holo.holoRequestNumber=ExpStruct.Holo.holoRequestNumber+1;
ExpStruct.Holo.holoRequests{ExpStruct.Holo.holoRequestNumber}=holoRequest;

% gen sequence list
Seq=makeHoloSequences(holoStimParams, setKey);
% make the daq sequences
makeHoloTrigSeqs(Seq, holoStimParams, holoRequest);
% save into exp struct
saveExpStructVars(holoRequest, holoStimParams);
 