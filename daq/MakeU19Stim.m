%% set params
startTime=1000; %ms; %typically 500
pulseDuration= 5; %ms Stimulation pulse 
TrigDuration = 5; %ms SLM flip command
stimFreq= 1; % Shouldn't matter but is used

P1=[0.055];%[0.01 0.025 0.05 0.1 .15]; %List of Powers To Use

% note: modded this to make all holos get the same number of pulses by
% default
stimTypeNumPulse = 20;
stimTypeHz = 30;

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
    powerList{c}       = [P1];
    waitList(c)        = 10; %time after stim before next stim ms
    hzList(c)          = stimTypeHz;  %30Hz
    pulseList(c)       = stimTypeNumPulse; %5 is default
    holosPerCycle(c)   = 1; %groups to interleave
    cellsPerHolo(c)    = numel(holosToUse{i});
    divTotalCells(c)   = numHolos*numel(stimTypeNumPulse); %divide total number of holos %i don't think it matters
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
%     totalCells = numel(unique([cellsToUse{:}]));
    temp2 = cellsToUse';
    totalCells = numel(unique(cat(1,temp2{:})));
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
holoStimParams.visStartTime = 1000;

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
    rois{i} =  cellsToUse{i};
    setKey{i} = i;
end

% for i = unique(holosets)
% 
% 
% these_rois=[];
% setlink_orders = cell(max(setlinks),1);
% for iset = unique(holoSets)
%     if length(unique(cellsPerHolo(holoSets==iset))) > 1 || length(unique(nHolos(holoSets==iset))) > 1 ||  length(unique(setlinks(holoSets==iset))) > 1
%         error('Holoset has multiple cellsPerHolo or nHolos');
%     else
%         nPerHolo = cellsPerHolo(find(holoSets==iset,1)); %TODO error if there are multiple cell numbers for sets
%         nHolo = nHolos(find(holoSets==iset,1));
%     end
%     %then create a unique + random order
%     this_set_order = 1:totalCells;%randperm(totalCells);
%     %check if this is part of a set
%     if setlinks(holoSets==iset)>0
%         this_setlink = setlinks(find(holoSets==iset,1));
%         %if already has an order, override the generated one
%         if ~isempty(setlink_orders{this_setlink})
%             this_set_order = setlink_orders{this_setlink};
%         end
%         %replace the setlink_orders with only the unused rois
%         setlink_orders{this_setlink} = this_set_order(nHolo*nPerHolo+1:end);
%     end
%     this_set_order = this_set_order(1:nHolo*nPerHolo);
% %     these_rois = randperm(totalCells,nPerHolo); %temp stand in
%     these_rois = makeHoloRois(nPerHolo,this_set_order);
%     setKey{iset} = [numel(rois)+1:numel(rois)+numel(these_rois)];
%     rois = [rois, cellfun(@(x) cellsToUse(x),these_rois,'uniformoutput',0)];
% end


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
 