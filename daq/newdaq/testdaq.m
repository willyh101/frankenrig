clear
close all
clc

addpath(genpath('.'))

fprintf('Starting daq...\r')

fprintf('Loading defaults... ')
setup = getDefaults();
pause(0.5)
fprintf('OK.\n')

fprintf('Making MATLAB NIDAQ object... ')
dq = daq('ni');
dq.Rate = setup.daqrate;
pause(0.5)
fprintf('OK.\n')

%% OUTPUTS

laser_eom = AnalogOutput('laser EOM', 2);
si_trig = DigitalOutput('si trig', 0, 0);
pt_trig = DigitalOutput('pt trig', 0, 3);
slm_trig = DigitalOutput('slm trig', 0, 2);

outputs = DAQoutput(dq, 'ouput');
outputs.addio(laser_eom, si_trig, pt_trig, slm_trig)
outputs.addToDaq()

%% INPUTS

slm_flip = AnalogInput('slm flip', 0);
% run_speed = EdgeCounter('run speed', )
run_speed = DigitalInput('run speed', 0, 10);
pt_stim_id = DigitalInput('pt stim id', 0, 13);
pt_stim_clk = DigitalInput('pt clock', 0, 14);
pt_flip = DigitalInput('pt flip', 0, 15);

inputs = DAQgroup(dq, 'input');
inputs.addio(slm_flip, run_speed, pt_stim_id, pt_stim_clk, pt_flip)
inputs.addToDaq()

%% create outputs
sweepLength = 2;

this_sweep = outputs.createSweep(sweepLength);

% ScanImage trigger
tg = outputs.createPulseOutput(sweepLength);
tg.addStartTrigger()
si_trig.output = tg;
this_sweep(:, si_trig) = tg;

% psych trigger
tg = outputs.createPulseOutput(sweepLength);
tg.addPulse(400)
pt_trig.output = tg;
this_sweep(:, pt_trig) = tg;

% slm trigger
slm_trig_times = [500:30:500+30*5, 700:30:700+30*5, 900:30:900+30*5];
tg = outputs.createPulseOutput(sweepLength);
tg.addPulseTrain(slm_trig_times)
slm_trig.output = tg;
this_sweep(:, slm_trig) = tg;

% laser EOM
laser_eom_times = [500:30:500+30*5, 700:30:700+30*5, 900:30:900+30*5]+300;
tg = outputs.createPulseOutput(sweepLength);
tg.addPulseTrain(laser_eom_times)
laser_eom.output = tg;
this_sweep(:, laser_eom) = tg;

figure(1)
clf
imagesc(this_sweep)


%% trial/condition manager

%% send/capture data

% to be turned into a function
indata = readwrite(dq, this_sweep, 'OutputFormat', 'Matrix');

%%





