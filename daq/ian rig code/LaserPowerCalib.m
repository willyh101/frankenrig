%function LaserPowerCalib


locations = FrankenScopeRigFile();

s = daq.createSession('ni'); %ni is company name

s.Rate=20000;
s.ExternalTriggerTimeout=30000000; %basically never time out
s.addAnalogOutputChannel('Dev1',2,'Voltage');

load(locations.PowerCalib,'LaserPower');


%%
disp('Do this calibration with the z-block in place')
disp('On the Holography Computer put a hologram near the zero order')
disp('Put Laser Gate on Bypass') %added 3/18/21

eomOffset = -0.15;
Fs = s.Rate;
trialLength = 0.10;%s

coarseCurve = linspace(eomOffset,5,25);
fullPowerCurve = coarseCurve(1:9);

extraFineCurve = linspace(eomOffset,coarseCurve(2),15); % was 15
fineCurvePart = linspace(extraFineCurve(end),coarseCurve(6),16);
fineCurve = [extraFineCurve fineCurvePart(2:end)];

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

    dataCoarse(i)=input(['Volt: ' num2str(Volt) ' Power (mW)? ']);
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

    dataFullPower(i)=input(['Volt: ' num2str(Volt) ' Power (mW) ? ']);
end

%Now doing low power
disp('Doing Low Power Curve')
dataLowPower =[];
for i=1:numel(fineCurve)
    
    Volt = fineCurve(i);
    
    if Volt==0
        Volt = 1e-5;
    end
    [output] = makepulseoutputs(1, 1, trialLength*Fs, Volt ,1, Fs, trialLength);
    output(output==0)=eomOffset;
    output(end-20:end)=[];
    s.queueOutputData([output]);
    s.startForeground;

    dataLowPower(i)=input(['Volt: ' num2str(Volt) ' Power (mW)? ']);
end

Volt=0;
[output] = makepulseoutputs(1, 1, trialLength*Fs, Volt ,1, Fs, trialLength);
output(output==0)=eomOffset;
output(end-20:end)=[];
s.queueOutputData([output]);
s.startForeground;

%%

convFactor = dataFullPower./dataCoarse(1:9);
convFactor = mean(convFactor(5:end));

dataEst = dataCoarse*convFactor;

figure(202);clf
plot(coarseCurve,dataEst,'-o');
hold on
plot(fineCurve,dataLowPower,'-o');

laserToSave = [fineCurve coarseCurve(7:end)];
dataToSave = ([dataLowPower dataEst(7:end)])./1000;

%% now you should enter the diffraction efficiency of your holo and divide by that
de = input('Now enter the diffraction efficiency of your holo!:');
dataToSave = dataToSave/de;
%% Ensure Monotonically increasing
while any(diff(dataToSave)==0)
    dataToSave(find(diff(dataToSave)==0)+1)=dataToSave(find(diff(dataToSave)==0,1)+1)+0.00001;
end
%% Overwrite old Laser Powers
LaserPower.EOMVoltage = laserToSave;
LaserPower.PowerOutputTF = dataToSave;
LaserPower.DEofTest = de;
LaserPower.dateCollected = date;
%%
save(locations.PowerCalib,'LaserPower');
%Save a Local copy
save(['C:\data\LaserPowerCalibs\LaserPower_' date '.mat'],'LaserPower');
save('C:\data\LaserPowerCalibs\LaserPower.mat','LaserPower');

disp('Saved')
figure(101)
plot(laserToSave,dataToSave)
disp(['Max power (raw) is: ' num2str(dataEst(end)/1000) ' Watts.'])
disp(['Max power (DE-adjusted) is: ' num2str(dataToSave(end)) ' Watts.'])
