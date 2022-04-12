% run a simple daq to test code coming in from a source
% going to use the new daq code from matlab
clear
clc

wait_time = 15;

dq = daq('ni');
dq.Rate = 10000;

% should be pt trig stim on/off
addinput(dq, 'Dev1', 0, 'Voltage')
% addinput(dq, 'Dev1', 'port0/line15', 'Digital'); % stim indicators
addinput(dq, 'Dev1', 'port0/line8:15', 'Digital'); % wide range here

disp('ready')
dataIn = read(dq, seconds(wait_time), 'OutputFormat', 'Matrix');
disp('done.')

figure(1)
clf
imagesc(dataIn)

%% compute the resulting stim info

dat = dataIn(:,6:7);
clk = dataIn(:,6);
stm = dataIn(:,7);
clk_times = find(diff(clk) == 1);
stm_times = find(diff(stm) == 1);
binvec = ismember(clk_times, stm_times)';
binaryVectorToDecimal(binvec)

%%
decodeBinaryDaqSignal(clk, stm)

%% test exporting clock to PT computer
clear
clc

dq = daq('ni');
dq.Rate = 10000;

addinput(dq, 'Dev1', 0, 'Voltage')
% addinput(dq, 'Dev1', 'port0/line15', 'Digital'); % stim indicators
% addinput(dq, 'Dev1', 'port0/line15', 'Digital'); % sample clock
addinput(dq, 'Dev1', 'port0/line8:15', 'Digital'); % all digital inputs








