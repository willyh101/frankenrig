
locations = FrankenScopeRigFile();

s = daq.createSession('ni'); %ni is company name

s.Rate=20000;
s.ExternalTriggerTimeout=30000000; %basically never time out
s.addAnalogOutputChannel('Dev1',2,'Voltage');

load(locations.PowerCalib,'LaserPower');

%% Simple request mW and measure power

PowerAsk = 100; % note this is in mW
DE = 0.623; % enter DE of requested hologram

eomOffset = -0.15;
Fs = s.Rate;
trialLength = 0.10;%s

PowerAsk = PowerAsk/1000;
PowerRequest = PowerAsk/DE;

Volt = function_EOMVoltage(LaserPower.EOMVoltage,LaserPower.PowerOutputTF,PowerRequest);
[output] = makepulseoutputs(1, 1, trialLength*Fs, Volt ,1, Fs, trialLength);
output(output==0)=eomOffset;
output(end-20:end)=[];

disp(['You are requesting ' num2str(PowerRequest) ' Watts of power.'])
disp(['The EOM voltage needed for this is ' num2str(Volt) ' Volts.'])
str = input('Are you sure? (y/n) ', 's');
if str ~= 'y'
    Volt=0;
    [output] = makepulseoutputs(1, 1, trialLength*Fs, Volt ,1, Fs, trialLength);
    output(output==0)=eomOffset;
    output(end-20:end)=[];
    s.queueOutputData([output]);
    s.startForeground;
    disp('...quitting.')
    return
else
    s.queueOutputData([output]);
    s.startForeground;
end

input('Press enter when done.')

Volt=0;
[output] = makepulseoutputs(1, 1, trialLength*Fs, Volt ,1, Fs, trialLength);
output(output==0)=eomOffset;
output(end-20:end)=[];
s.queueOutputData([output]);
s.startForeground;
%%
disp('Do this calibration with the z-block in place')
disp('On the Holography Computer put a hologram near the zero order')

eomOffset = -0.15;
Fs = s.Rate;
trialLength = 0.10;%s

coarseCurve = linspace(eomOffset,5,20);
fullPowerCurve = coarseCurve(1:5);

disp('First a Full Power Curve in divided Mode')
disp('Use divided Mode 10. With the Laser Gate Turned Off')
dataCoarse =[];
for i=1:numel(coarseCurve)
    
    Volt = coarseCurve(i);
    
    if Volt==0
        Volt = 1e-5;
    end
    [output] = makepulseoutputs(1, 1, trialLength*Fs, Volt ,1, Fs, trialLength);
    output(output==0)=eomOffset;
    output(end-20:end)=[];
    s.queueOutputData([output]);
    s.startForeground;

    dataCoarse(i)=input(['Volt: ' num2str(Volt) ' Power? ']);
end

Volt=0;
[output] = makepulseoutputs(1, 1, trialLength*Fs, Volt ,1, Fs, trialLength);
output(output==0)=eomOffset;
output(end-20:end)=[];
s.queueOutputData([output]);
s.startForeground;

dummy = input('Now put the Laser on nondivided full power mode (press enter when ready)');
dataFullPower =[];
for i=1:numel(fullPowerCurve)
    
    Volt = fullPowerCurve(i);
    
    if Volt==0
        Volt = 1e-5;
    end
    [output] = makepulseoutputs(1, 1, trialLength*Fs, Volt ,1, Fs, trialLength);
    output(output==0)=eomOffset;
    output(end-20:end)=[];
    s.queueOutputData([output]);
    s.startForeground;

    dataFullPower(i)=input(['Volt: ' num2str(Volt) ' Power? ']);
end

Volt=0;
[output] = makepulseoutputs(1, 1, trialLength*Fs, Volt ,1, Fs, trialLength);
output(output==0)=eomOffset;
output(end-20:end)=[];
s.queueOutputData([output]);
s.startForeground;

%%

convFactor = dataFullPower./dataCoarse(1:5);
convFactor = mean(convFactor(5:end));

dataEst = dataCoarse*convFactor;

figure(202);clf
plot(coarseCurve,dataEst,'-o');

laserToSave = coarseCurve;
dataToSave = dataEst./1000;

%% now you should enter the diffraction efficiency of your holo and divide by that
de = input('Now enter the diffraction efficiency of your holo!:');
dataToSave = dataToSave/de;
%% Ensure Monotonically increasing
while any(diff(dataToSave)==0)
    dataToSave(find(diff(dataToSave)==0)+1)=dataToSave(find(diff(dataToSave)==0,1)+1)+0.00001;
end
%%
figure(101)
plot(laserToSave,dataToSave)


%% JUST GET MAX POWER
disp('About to find max power. QUADRUPLE CHECK THE LASER IS IN 10 DIV MODE or you will have to buy a new SLM!')
input('When you are sure, press ENTER')
str = input('Last chance to abort (q)...', 's');
if str == 'q'
    return
else
    eomOffset = -0.15;
    Fs = s.Rate;
    trialLength = 0.10;%s
    Volt = 1;
    [output] = makepulseoutputs(1, 1, trialLength*Fs, Volt ,1, Fs, trialLength);
    output(output==0)=eomOffset;
    output(end-20:end)=[];
    
    s.queueOutputData([output]);
    s.startForeground;
    str = input('this is 1 V, does this make sense for div mode 10 (y/n): ', 's');
    if str ~= 'y'
        Volt=0;
        [output] = makepulseoutputs(1, 1, trialLength*Fs, Volt ,1, Fs, trialLength);
        output(output==0)=eomOffset;
        output(end-20:end)=[];
        s.queueOutputData([output]);
        s.startForeground;
        return
    else
        eomOffset = -0.15;
        Fs = s.Rate;
        trialLength = 0.10;%s
        Volt = 3;
        [output] = makepulseoutputs(1, 1, trialLength*Fs, Volt ,1, Fs, trialLength);
        output(output==0)=eomOffset;
        output(end-20:end)=[];
        s.queueOutputData([output]);
        s.startForeground;
        pwr = input('Power: ');
        
        Volt=0;
        [output] = makepulseoutputs(1, 1, trialLength*Fs, Volt ,1, Fs, trialLength);
        output(output==0)=eomOffset;
        output(end-20:end)=[];
        s.queueOutputData([output]);
        s.startForeground;
        
        de = input('enter DE of hologram: ');
        maxpwr = pwr/de;
        disp(['max power (DE adjusted) is:  ' num2str(maxpwr)])
    end
end

