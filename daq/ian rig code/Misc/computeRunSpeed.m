function [running_speed] = computeRunSpeed(MotorA,binWidth)

% default rotary encoder has 360 pulses per revolution
% default circumference is 2*pi*7.6 cm = 47.75cm/revolution or 18.84 inches/revolution

% global Exp_Defaults

Fs = 20000; %hardcoded frame speed, if you ever change it it will stop working

%speed is determined by MotorA right now
runTicks = [double(diff(MotorA)>0); 0];
% tickTimes = find(MotorA==0);

%binWidth = 0.1; % in s

binInPnts = round(Fs * binWidth) ; %set bin size to frame rate

binnedTicks=[];
nBins = length(runTicks) /binInPnts;
for i=1:nBins
    binnedTicks(i)  = sum(runTicks((i-1)*binInPnts+1:i*binInPnts));
end


circumference = 47.75; %cm at edge of wheel
angularDistance = binnedTicks/360;
LinearDistance = angularDistance * circumference;

running_speed = LinearDistance / binWidth;
