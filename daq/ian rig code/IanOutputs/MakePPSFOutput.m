%% if desired, load hr
locations = FrankenScopeRigFile();
load([locations.HoloRequest_DAQ 'holoRequest.mat']);
disp('loaded holorequest!')
disp(['This holorequest has ' num2str(size(holoRequest.targets,1)) ' targets.'])
%% Determine Stimmable
out = ExpStruct.stimTest;

PSTHs=out.PSTHs;
powers = out.powers;
holoRequest.onlinePowerCurve = out;
refPower = 50/1000;
[hr]=roiWeightsFromPSTHs(PSTHs,powers,holoRequest,refPower);
% eligibleCells = find(~isnan(hr.roiWeights)); %find(hr.roiWeights==1);
noResponseCells = find(isnan(hr.roiWeights));
eligibleCells = find(hr.roiWeights<=1.5);

ExpStruct.eligibleCells = eligibleCells;
ExpStruct.noResponseCells = noResponseCells;
numel(eligibleCells)
holoRequest.roiWeights = hr.roiWeights;

%% Select inital Targets
uniqueDepths = unique(holoRequest.targets(:,3));

eligibleCells(holoRequest.targets(eligibleCells,3)>uniqueDepths(2))=[];
numel(eligibleCells)

numCellsPerHolo = 10; 
cellsToUse = spreadOutTargetCells(eligibleCells,numCellsPerHolo,holoRequest,0); 
holosToUse= {cellsToUse(1:numCellsPerHolo)};

holosToUse{2} = cellsToUse(numCellsPerHolo+1:numCellsPerHolo*2);

%% Set PPSF Incriments
pxPerMu = 512/800; %pixels per um 
optPerMu = 55/60; %optotune units per um

Xvals = (-3:3:30) * pxPerMu; %offsets in um
Zvals = (-6:6:60) * optPerMu; %offsets in um

numTargets = size(holoRequest.targets,1);

targetsOrig = holoRequest.targets;
roiWeightsOrig = holoRequest.roiWeights;
holosToUseOrig = holosToUse; 

holoAddedOffsets=[];
holosToUse=[];
for n = 1:numel(holosToUseOrig)
    for i = 1:numel(Xvals)
        x =  targetsOrig +  [Xvals(i) 0 0];
        thisHoloSet = holosToUseOrig{n};
        thisHoloSet(x(thisHoloSet,1)<1 | x(thisHoloSet,1)>511)=[];
        
        holoRequest.targets = cat(1, holoRequest.targets,x);
        holoRequest.roiWeights = cat(1, holoRequest.roiWeights,roiWeightsOrig);
        
        holosToUse{end+1} = thisHoloSet+numTargets*i;
        
        holoAddedOffsets{end+1} = [Xvals(i) 0 0];
    end
    
    for i = 1:numel(Zvals)
        x =  targetsOrig +  [0 0 Zvals(i)];
        thisHoloSet = holosToUseOrig{n};
        thisHoloSet(x(thisHoloSet,3)<0 | x(thisHoloSet,3)>80)=[];
        
        holoRequest.targets = cat(1, holoRequest.targets,x);
        holoRequest.roiWeights = cat(1, holoRequest.roiWeights,roiWeightsOrig);
        holosToUse{end+1} = thisHoloSet+numTargets*i;
        
        holoAddedOffsets{end+1} = [0 0 Zvals(i)];
        
    end
end
holoStimParams.holoAddedOffsets = holoAddedOffsets; 

%% set params
startTime=1000; %ms; %typically 500
pulseDuration= 5; %ms Stimulation pulse 
TrigDuration = 5; %ms SLM flip command
stimFreq= 1; % Shouldn't matter but is used

P1=0.055 ; %[0.055];%[0.01 0.025 0.05 0.1 .15]; %List of Powers To Use

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
    waitList(c)        = 10; %time after stim before next stim ms
    hzList(c)          = 10; %stimTypeHz(k);  %30Hz
    pulseList(c)       = 10; %stimTypeNumPulse(k); %5 is default
    holosPerCycle(c)   = 1; %groups to interleave
    cellsPerHolo(c)    = numel(holosToUse{i});
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
    try ;
        totalCells = numel(unique(cat(1,cellsToUse{:})));
    catch
        totalCells = numel(unique(cat(2,cellsToUse{:})));
    end
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

setKey = {};
rois = {};
for i = 1:numHolos%*numel(stimTypeNumPulse)
    rois{i} =  cellsToUse{i};
    setKey{i} = i;
end

holoRequest.rois = rois;

%% run the universal stuff
% talk to holo comp, fill in roi weights if its not there
msocketPrep;
holoRequest = transferHR(holoRequest);

ExpStruct.Holo.holoRequestNumber=ExpStruct.Holo.holoRequestNumber+1;
ExpStruct.Holo.holoRequests{ExpStruct.Holo.holoRequestNumber}=holoRequest;
%%E
% gen sequence list-
Seq=makeHoloSequences(holoStimParams, setKey); 
% make the daq sequences
makeHoloTrigSeqs(Seq, holoStimParams, holoRequest);
% save into exp struct
saveExpStructVars(holoRequest, holoStimParams);
 