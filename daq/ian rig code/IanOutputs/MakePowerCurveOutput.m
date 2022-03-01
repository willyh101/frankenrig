%% if desired, load hr
locations = FrankenScopeRigFile();
load([locations.HoloRequest_DAQ 'holoRequest.mat']);
disp('loaded holorequest!')
disp(['This holorequest has ' num2str(size(holoRequest.targets,1)) ' targets.'])
%% set params
startTime=1000; %ms; %set to 1s on 2/27/20 (previously 500ms) -Ian
pulseDuration = 5; %ms Stimulation pulse 
TrigDuration = 5; %ms SLM flip command
stimFreq= 1; % Shouldn't matter but is used

P1= [0.0125 0.025 0.05]; %[.025 .05 .075 .1]; %[.01 .02 .03 .04 .050 .075];%[.05 .06 .07 .08 .09 .1];%[.01 .025 .04 .055 .07 .085 0.1 0.115 0.13]; %[0.05 0.1 0.125 0.15];%[.01 .02 .03 .04 .05];%[0.05 .1 .15];%List of Powers To Use
cellsToUse = 0;%ExpStruct.eligibleCells; find(holoRequest.targets(:,3)<75); %0;%[12 13 6 15 16 14 ]; %[6 12 13 14 15 16]; %0; %Put List of Targets to shoot here, or zero to use all available.


if any(P1>0.3)
    errordlg('Did you want >300mW???')
end

powerList       = { P1 };
waitList        = [ 6 ]; %time after stim before next stim ms
hzList          = [ 30 ];
pulseList       = [ 5 ]; %previously 10; 5 is standard as of 5/15/19
holosPerCycle   = [ 1 ]; %groups to interleave
cellsPerHolo    = [ 1 ]; %number of cells per holo
divTotalCells   = [ 1 ]; %divide total number of holos 
holoSets        = [ 1 ]; %unique groups of cells per 
setlinks        = [ 1 ];

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
% gen sequence list-
Seq=makeHoloSequences(holoStimParams, setKey); 
% make the daq sequences
makeHoloTrigSeqs(Seq, holoStimParams, holoRequest);
% save into exp struct
saveExpStructVars(holoRequest, holoStimParams);
 