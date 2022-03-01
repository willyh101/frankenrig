function [runVector] = runReader2(ExpStruct,epoch,FR,chan)
if nargin<4
    chan = 5;
end

sws=[];
sws(1) = ExpStruct.EpochEnterSweep{epoch};
try
    sws(2) = ExpStruct.EpochEnterSweep{epoch+1}-1;
catch
    sws(2) = ExpStruct.sweep_counter-1;
end
sws = sws(1):sws(2);


digitSwp = ExpStruct.digitalSweeps(sws);%% this is where run speed is
digitSwp = cellfun(@(x) x(:,chan),digitSwp,'uniformoutput',0); %select line 1 

% digitSwp = cellfun(@(x) x(:,1),digitSwp,'uniformoutput',0); %select line 1 

runVector = cellfun(@(x) computeRunSpeed(x,1/FR),digitSwp,'uniformoutput',0); %convert to cm/s

runVector=cell2mat(runVector');


