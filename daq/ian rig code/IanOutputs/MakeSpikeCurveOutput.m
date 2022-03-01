%% if desired, load hr
locations = FrankenScopeRigFile();
load([locations.HoloRequest_DAQ 'holoRequest.mat']);
disp('loaded holorequest!')
%% set params
startTime=1000; %ms;
pulseDuration = 5; %ms Stimulation pulse 
TrigDuration = 5; %ms SLM flip command
stimFreq= 1; % Shouldn't matter but is used

numCellsPerHolo = 5;% Set this to change number of simultaneous cells

P1= [.06];%[0.05 .1 .15];%List of Powers To Use
S1 = [1 3 5 7 10 15 20]; %updated 7/28/2020
% cellsToUse = spreadOutTargetCells(look_right,numCellsPerHolo,holoRequest,0);%look_right; %cellsToCalib;%   0;%[1:10]; %Put List of Targets to shoot here, or zero to use all available.
cellsToUse = spreadOutTargetCells(cellsToCalib,numCellsPerHolo,holoRequest,0);%look_right; %cellsToCalib;%   0;%[1:10]; %Put List of Targets to shoot here, or zero to use all available.

if any(P1>0.3)
    errordlg('Did you want >300mW???')
end


powerList       = {  };
waitList        = [ ]; %time after stim before next stim ms
hzList          = [  ];
pulseList       = [  ]; %previously 10
holosPerCycle   = [  ]; %groups to interleave
cellsPerHolo    = [  ];
divTotalCells   = [  ]; %divide total number of holos 
holoSets        = [  ]; %unique groups of cells per 
setlinks        = [  ];
bwnGroupPause   = [ 10 ]; %pause between groups in ms (default 10ms for normal 250 for detailed


for i=1:numel(S1)
    powerList(i)       = { P1 };
    waitList(i)        = [ 6 ]; %time after stim before next stim ms
    hzList(i)          = S1(i);
    pulseList(i)       = S1(i); %previously 10; 5 is standard as of 5/15/19
    holosPerCycle(i)   = [ 1 ]; %groups to interleave
    cellsPerHolo(i)    = [ numCellsPerHolo ];  %num Cells ata  time
    divTotalCells(i)   = [ 1 ]; %divide total number of holos
    holoSets(i)        = [ 1 ]; %unique groups of cells per
    setlinks(i)        = [ 1 ];
end
% for i=1:3
%     powerList{i}       = [P1];
%     waitList(i)        = 10; %time after stim before next stim ms
%     hzList(i)          = 50;
%     pulseList(i)       = 10;
%     holosPerCycle(i)   = 1; %groups to interleave
%     cellsPerHolo(i)    = 1;
%     divTotalCells(i)   = 3; %divide total number of holos
%     holoSets(i)        = i; %unique groups of cells per
%     setlinks(i)        = i;
% end
if cellsToUse==0
totalCells =size(holoRequest.targets,1);
disp(['Total Cells Detected ' num2str(totalCells)]);
cellsToUse = 1:totalCells;
else
    totalCells = numel(cellsToUse);
    disp(['Using ' num2str(totalCells) ' Cells']);
end

nHolos          = floor(totalCells./divTotalCells./cellsPerHolo); %only make complete holograms
repsList        = floor(nHolos./holosPerCycle); 

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
holoStimParams.nHolos = nHolos;
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
rois = {};
setKey = {};
these_rois=[];
setlink_orders = cell(max(setlinks),1);
for iset = unique(holoSets)
    if length(unique(cellsPerHolo(holoSets==iset))) > 1 || length(unique(nHolos(holoSets==iset))) > 1 ||  length(unique(setlinks(holoSets==iset))) > 1
        error('Holoset has multiple cellsPerHolo or nHolos');
    else
        nPerHolo = cellsPerHolo(find(holoSets==iset,1)); %TODO error if there are multiple cell numbers for sets
        nHolo = nHolos(find(holoSets==iset,1));
    end
    %then create a unique + random order
    this_set_order = 1:totalCells;%randperm(totalCells);
    %check if this is part of a set
    if setlinks(holoSets==iset)>0
        this_setlink = setlinks(find(holoSets==iset,1));
        %if already has an order, override the generated one
        if ~isempty(setlink_orders{this_setlink})
            this_set_order = setlink_orders{this_setlink};
        end
        %replace the setlink_orders with only the unused rois
        setlink_orders{this_setlink} = this_set_order(nHolo*nPerHolo+1:end);
    end
    this_set_order = this_set_order(1:nHolo*nPerHolo);
%     these_rois = randperm(totalCells,nPerHolo); %temp stand in
    these_rois = makeHoloRois(nPerHolo,this_set_order);
    setKey{iset} = [numel(rois)+1:numel(rois)+numel(these_rois)];
    rois = [rois, cellfun(@(x) cellsToUse(x),these_rois,'uniformoutput',0)];
end
holoRequest.rois = rois;

%% run the universal stuff
% talk to holo comp, fill in roi weights if its not there
msocketPrep;
holoRequest = transferHR(holoRequest);

ExpStruct.Holo.holoRequestNumber=ExpStruct.Holo.holoRequestNumber+1;
ExpStruct.Holo.holoRequests{ExpStruct.Holo.holoRequestNumber}=holoRequest;
%%E
% gen sequence list
Seq=makeHoloSequences(holoStimParams, setKey); 
% make the daq sequences
makeHoloTrigSeqs(Seq, holoStimParams, holoRequest);
% save into exp struct
saveExpStructVars(holoRequest, holoStimParams);
 