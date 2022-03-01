function testStroboscopicDAQ

locations = FrankenScopeRigFile();

s = daq.createSession('ni'); %ni is company name

s.Rate=20000;
s.ExternalTriggerTimeout=30000000; %basically never time out
s.addAnalogOutputChannel('Dev3',2,'Voltage'); %Laser EOM
addDigitalChannel(s,'Dev3','Port0/Line2','OutputOnly');  %change Holos Line


load(locations.PowerCalib,'LaserPower');


%initalize contact
IP='128.32.173.87';

[HoloSocket]=msocketPrep(IP);

%%
%laser params
interval =0.5; %in ms
duration = 0.5; %in ms

trialLength = 0.075;%s
eomOffset = -0.15;
PowerRequest=0.3; %in watts %don't use divided mode 300mW about right

numTrials = trialLength*1000/interval;

numHolos = 4;
holoRate = 60; %in Hz;
holoDur =3;

OutputList=[];

for i=1:numTrials
    
    
    
    Volt = function_EOMVoltage(LaserPower.EOMVoltage,LaserPower.PowerOutputTF,PowerRequest);
    Fs = s.Rate;
    
    
    [output] = makepulseoutputs(1+(i-1)*interval, 1, duration, Volt ,1, Fs, trialLength);
    output(output==0)=eomOffset;
    laserOutput = output;
    
    [output] = makepulseoutputs(10, numHolos, holoDur, 1 ,holoRate, Fs, trialLength);
    output(output==0)=eomOffset;
    OutputList{i}= [laserOutput output ];
end

%%
nReps = 5;

for k=1:nReps
    for i=1:numTrials
        invar =[];
        while isempty(invar)
            invar = msrecv(HoloSocket,0.5);
        end
        if ~strcmp(invar,'ready')
            disp('Unexpected input stopping...')
            break
        else
            invar='flush';
            while ~isempty(invar)
                invar  = msrecv(HoloSocket,0.5);
            end
        end
        
        disp(['Cycle ' num2str(i) ' of ' num2str(numTrials) ' Rep ' num2str(k)])
        order = 1:numHolos;
        mssend(HoloSocket,order);
        
        s.queueOutputData(OutputList{i});
        s.startForeground;
        mssend(HoloSocket,[i k]);
    end
    
end
pause(1);
mssend(HoloSocket,'end');
disp('end')
