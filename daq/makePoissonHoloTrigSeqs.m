function makePoissonHoloTrigSeqs(Seq,holoStimParams,holoRequest)
global ExpStruct Exp_Defaults
locations = FrankenScopeRigFile();
load(locations.PowerCalib,'LaserPower');

disp('Making Sequences...')
DE_list = holoRequest.DE_list;

% ScanImage Trigger
ExpStruct.triggerSI5=zeros(size(ExpStruct.triggerSI5))...
    + makepulseoutputs(1,1,25,1,1,Exp_Defaults.Fs,size(ExpStruct.StimLaserEOM,1)/Exp_Defaults.Fs);

updateAOaxes();

blankOutput=zeros(size(ExpStruct.StimLaserEOM));

%independent vis start time control added 4/5/21 -IAO
% moved up here to fix the blank condition
if ~isfield(holoStimParams,'visStartTime')
    holoStimParams.visStartTime = holoStimParams.startTime-80;
end

%Make A blank output
ExpStruct.triggerPuffer=blankOutput;%ExpStruct.nextsequenceTrigger=blankOutput;
ExpStruct.nextholoTrigger = blankOutput;
ExpStruct.StimLaserEOM=blankOutput;
ExpStruct.nextsequenceTrigger=zeros(size(ExpStruct.nextsequenceTrigger))...
    + makepulseoutputs(holoStimParams.visStartTime,1,10,1,50,Exp_Defaults.Fs,size(ExpStruct.StimLaserEOM,1)/Exp_Defaults.Fs);

updateAOaxes();
saveoutput('Out 0. 0mW. No Pulses');
disp('Saved Control Output');

ExpStruct.bigListOfFirstStimTimes=[];



% p=1;
c = 0; %count the output number
newSeq=[];
for i = 1:numel(Seq) %Each type of stim
    t=tic;
    for p = 1:numel(holoStimParams.powerList{i})
        c=c+1;
        
        %Specify seq list
        ExpStruct.triggerPuffer=zeros(size(ExpStruct.triggerPuffer))...
            + makepulseoutputs(25,c,1.5,1,200,Exp_Defaults.Fs,size(ExpStruct.StimLaserEOM,1)/Exp_Defaults.Fs);
        %changed to 1.5ms at 200Hz on 4/20/19
        
        %Trig Vis
        ExpStruct.nextsequenceTrigger=zeros(size(ExpStruct.nextsequenceTrigger))...
            + makepulseoutputs(holoStimParams.visStartTime,1,10,1,50,Exp_Defaults.Fs,size(ExpStruct.StimLaserEOM,1)/Exp_Defaults.Fs);
        
        targetList=(Seq{i});
        
        for k=1:numel(targetList)
            targets=holoRequest.rois{targetList(k)};
            %             PowerRequest = (powerList{i}(p)*numel(targets))/DE_list(Seq{i}(k));
            
            %fix nan weights
            theseWeights = holoRequest.roiWeights(targets);
            theseWeights(isnan(theseWeights)) = 1;
            powerAsk = holoStimParams.powerList{i}(p)*sum(theseWeights);
            if DE_list(Seq{i}(k))==0
                disp('DE is 0. setting to 0V. probably no holo')
                Volt=0;
            else
                PowerRequest = powerAsk/DE_list(Seq{i}(k));
                
                
                Volt = function_EOMVoltage(LaserPower.EOMVoltage,LaserPower.PowerOutputTF,PowerRequest);
                if isnan(Volt)
                    disp('Could not set voltage picked 5 Volts')
                    Volt =5;% function_EOMVoltage(LaserPower.EOMVoltage,LaserPower.PowerOutputTF,max(LaserPower.PowerOutputTF));
                end
            end
            voltList(k)=Volt;
        end
        stimOutput = blankOutput;
        trigOutput = blankOutput;
        
        tm = holoStimParams.startTime;
        pulseStart = holoStimParams.startTime - holoStimParams.TrigDuration;
        counter = 1;
        Hz = holoStimParams.hzList(i);
        for R = 1:holoStimParams.repsList(i) %Repeat with a different number of cells
            try
                pulseStart = pulseStart + holoStimParams.bwnGroupPause;
            catch
            end
            
            tt = holoStimParams.poissonStimLength;
            minsep = holoStimParams.poissonMinSep;
            spks = poissonSpikeGenerator(Hz,tt,minsep);
            nspks = length(spks);
            
            for Pulse = 1:nspks
                tm = pulseStart + spks(Pulse);
                
                for Ce = 1:holoStimParams.holosPerCycle(i)
                    V = voltList(counter);
                    
                    % Trig Time
                    TT=makepulseoutputs(tm,1,...
                        holoStimParams.TrigDuration,...
                        1,holoStimParams.stimFreq,...
                        Exp_Defaults.Fs,...
                        size(ExpStruct.StimLaserEOM,1)/Exp_Defaults.Fs);
                    tm=tm+holoStimParams.waitList(i)+ ...
                        holoStimParams.TrigDuration;
                    trigOutput = trigOutput +TT;
                    
                    
                    %Stim Time
                    ST=makepulseoutputs(tm,1,holoStimParams.pulseDuration,V,holoStimParams.stimFreq,Exp_Defaults.Fs,size(ExpStruct.StimLaserEOM,1)/Exp_Defaults.Fs);
                    tm=tm+holoStimParams.pulseDuration-holoStimParams.TrigDuration;
                    stimOutput = stimOutput +ST;
                    

                end
            end
            counter=counter+1;
        end
        
        %if getFirstStimTimes breaks check that sweep is long enough
        try
            ExpStruct.bigListOfFirstStimTimes(:,c) = getFirstStimTimes(Seq{i}, trigOutput, holoRequest, Exp_Defaults.Fs, holoStimParams.totalCells, holoStimParams.pulseList(i));
        catch
            
            disp('WARNING: Unable to create a bigListOfFirstStimTimes')
            ExpStruct.bigListOfFirstStimTimes(:,c) = 0;
        end
        
        ExpStruct.nextholoTrigger = trigOutput;
        ExpStruct.StimLaserEOM=stimOutput;
        updateAOaxes();
        
        %         newSeq{c}=Seq{i};
        ExpStruct.bigListofSequences{c} = Seq{i};
        ExpStruct.powerList(c) = holoStimParams.powerList{i}(p);
        ExpStruct.pulseList(c) = holoStimParams.pulseList(i);
        ExpStruct.hzList(c)    = holoStimParams.hzList(i);
        % %         sendThis = ExpStruct.bigListofSequences{seqNum}
        % % sendThisSI.times = ExpStruct.bigListOfFirstStimTimes(:,seqNum);
        % % sendThisSI.power = ExpStruct.powerList(seqNum);
        
        saveoutput(['Out ' num2str(c) '. ' ...
            num2str(holoStimParams.powerList{i}(p)*1000) ...
            'mW ' num2str(Hz) 'Hz x' ...
            num2str(holoStimParams.pulseList(i)) ' ' ...
            num2str(holoStimParams.waitList(i)) 'ms wt ' ...
            num2str(holoStimParams.cellsPerHolo(i)) ' cell']);
        disp(['Saved power ' num2str(p) ' of ' ...
            num2str(numel(holoStimParams.powerList{i}))]);
    end
    disp(['Sequence ' num2str(i) ' complete. Took ' num2str(toc(t)) 's']);
end
