
%% set params
startTime=1000; %ms; %typically 500
pulseDuration= 5; %ms Stimulation pulse 
TrigDuration = 0.1; %ms SLM flip command
stimFreq= 1; % Shouldn't matter but is used

P1=[0.055];%[0.01 0.025 0.05 0.1 .15]; %List of Powers To Use

stimTypeNumTarget = unique(cellfun(@(x) numel(x),holosToUse)); %[3 10 25];%;[3 10];% 33];%[4 10];%[10 4 25];% 11];%[10 5 20]; %also dictates grouping needs to not be redundant
% stimTypeNumPulse = round(100./stimTypeNumTarget);%[33 20 10 5 3];%;[33 10];% 3]; %[25 10];% [10 25 4];% 10];%[10 20 5];
% stimTypeHz = stimTypeNumPulse;%[10 33 3];%10; [30 9];% 2.7];%[30 12]; %[12 30 4.8];% 12];%[15 30 7.5]; %added 10/3/19 to make all stims same length

numHolos=numel(holosToUse);

powerList       = {  };
waitList        = [ ]; %time after stim before next stim ms
hzList          = [  ];
pulseList       = [  ]; %previously 10
holosPerCycle   = [  ]; %groups to interleave
cellsPerHolo    = [  ];
divTotalCells   = [  ]; %divide total number of holos 
holoSets        = [  ]; %unique groups of cells per 
setlinks        = [  ];
c=1;

% for k= 1:numel(stimTypeNumPulse)
for i=1:numHolos
    k = find(stimTypeNumTarget==numel(holosToUse{i}));
    powerList{c}       = [P1];
    waitList(c)        = 5; %time after stim before next stim ms
    hzList(c)          = 10; %stimTypeHz(k);  %30Hz
    pulseList(c)       = 10; %stimTypeNumPulse(k); %Number of Pulses
    holosPerCycle(c)   = numel(holosToUse{i});%1; %groups to interleave
    cellsPerHolo(c)    = 1; %numel(holosToUse{i});
    divTotalCells(c)   = numHolos*numel(stimTypeNumTarget);%numHolos*numel(stimTypeNumPulse); %divide total number of holos %i don't think it matters
    holoSets(c)        = c; %unique groups of cells per
    setlinks(c)        = 1;
    c=c+1;
end
% end
cellsToUse =holosToUse;%repmat(holosToUse,[1 numel(stimTypeNumPulse)]);%[holosToUse{:}];% holosToUse;%GroupList;%cellsToShoot;eligibleCells;% 0; %Put List of Targets to shoot here, or zero to use all available.
% cellsToUse =[cellsToUse{:}];%
% cellsToUse =[holosToUse{:}];

if ~iscell(cellsToUse) && cellsToUse==0
    totalCells =size(holoRequest.targets,1);
    disp(['Total Cells Detected ' num2str(totalCells)]);
    cellsToUse = 1:totalCells;
elseif iscell(cellsToUse)
    totalCells = numel(unique(cat(1,cellsToUse{:})));
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
holoStimParams.divTotalCells = divTotalCells; %not used in subsequent functions
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

holoRequest.holoStimParams = holoStimParams;
%%gen holo rois


% if size(cellsToUse,1) > max(holoSets) 
%     disp('***ERROR*** not enough hologram targets')
% end
% 
% if size(cellsToUse,2) > max(cellsPerHolo)
%     disp('***ERROR*** not enough cells per hologram')
% end

setKey = {};
rois = {};
for i = 1:numHolos%*numel(stimTypeNumPulse)
    setKey{i} = cellsToUse{i}; %setKey points each pulse to an ROI
end

allCellsUnrolled = [cellsToUse{:}];
for i = 1:max(allCellsUnrolled(:))
    rois{i}= i; %code each ROI to match a target
end

holoRequest.rois = rois;

%% run the universal stuff
% talk to holo comp, fill in roi weights if its not there
msocketPrep;
holoRequest = transferHR(holoRequest);

ExpStruct.Holo.holoRequestNumber=ExpStruct.Holo.holoRequestNumber+1;
ExpStruct.Holo.holoRequests{ExpStruct.Holo.holoRequestNumber}=holoRequest;

% gen sequence lista
Seq=makeHoloSequences(holoStimParams, setKey);
% make the daq sequences
makeHoloTrigSeqs(Seq, holoStimParams, holoRequest);
% save into exp struct
saveExpStructVars(holoRequest, holoStimParams);
 