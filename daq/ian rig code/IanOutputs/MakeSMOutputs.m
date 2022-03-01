%% if desired, load hr
locations = FrankenScopeRigFile();
load([locations.HoloRequest_DAQ 'holoRequest.mat']);


%% set params
startTime=1000; %ms;
pulseDuration=5; %ms Stimulation pulse 
TrigDuration = 5; %ms SLM flip command
stimFreq= 1; % Shouldn't matter but is used

P1=[0.06]; %[0.01 0.025 0.05 .075 0.1 .15]; %List of Powers To Use

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
    pulseList(i)       = 10; %5 is default
    holosPerCycle(i)   = 1; %groups to interleave
    cellsPerHolo(i)    = numel(holosToUse{i});
    divTotalCells(i)   = numHolos; %divide total number of holos
    holoSets(i)        = i; %unique groups of cells per
    setlinks(i)        = 1;
end
cellsToUse = holosToUse;%GroupList;%cellsToShoot;eligibleCells;% 0; %Put List of Targets to shoot here, or zero to use all available.
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

holoRequest.holoStimParams = holoStimParams;

setKey = {};
rois = {};
for i = 1:numHolos
    rois{i} = cellsToUse{i}; %cellsToUse(i,:);
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
%%
% gen sequence list
Seq=makeHoloSequences(holoStimParams, setKey);
% make the daq sequences
makeHoloTrigSeqs(Seq, holoStimParams, holoRequest);
% save into exp struct
saveExpStructVars(holoRequest, holoStimParams);
 