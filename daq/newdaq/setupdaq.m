clear
close all
clc

fprintf('Starting daq...\r')
exp_datetime = datetime('now');
daqStartDate = datestr(exp_datetime, 'mm/dd/yyyy');
daqStartTime = datestr(exp_datetime, 'HH:MM:SS');
pause(0.5)

fprintf('Loading defaults... ')
setup = getDefaults();
pause(0.5)
fprintf('OK.\n')

fprintf('Making MATLAB NIDAQ object... ')
dq = daq('ni');
dq.Rate = setup.daqrate;
pause(0.5)
fprintf('OK.\n')

fprintf('Setting up IO channels for NIDAQ...')

% Inputs
inputs.names = {'slm flip', 'running', 'pt stim id', 'pt stim clk', 'pt flip'};
inputs.locs = {'ai0', 'port0/line10', 'port0/line13', 'port0/line14', 'port0/line15'};
assert(numel(inputs.names)==numel(input.locs), 'Inconsistent number of input names and channel locations.');

addinput(dq, setup.device, 0, 'Voltage');  % slm flip
addinput(dq, setup.device, 'port0/line10', 'Digital'); % running
addinput(dq, setup.device, 'port0/line13', 'Digital'); % pt stim id
addinput(dq, setup.device, 'port0/line14', 'Digital'); % pt stim id clock signal
addinput(dq, setup.device, 'port0/line15', 'Digital'); % pt stim indicator

% Outputs
outputs.names = {'laser EOM', 'si trig', 'slm trig', 'pt trig'};
outputs.locs = {'ao2', 'port0/line0', 'port0/line2', 'port0/line3'};
assert(numel(outputs.names)==numel(outputs.locs), 'Inconsistent number of output names and channel locations.');

addoutput(dq, setup.device, 2, 'Voltage'); % laser EOM
addoutput(dq, setup.device, 'port0/line0', 'Digital'); % si trig
addoutput(dq, setup.device, 'port0/line2', 'Digital'); % slm trig
addoutput(dq, setup.device, 'port0/line3', 'Digital'); % pt trig

pause(0.5)
fprinf('OK.\n')

fprintf('\n#---Inputs Added---#\n')
for i=1:numel(inputs.names)
    fprintf([num2str(i) ': ' inputs.names{i} '\n']);
    pause(0.1)
end
pause(1)

fprintf('\n#---Outputs Added---#\n')
for i=1:numel(inputs)
    fprintf([num2str(i) ': ' outputs.names{i} '\n']);
    pause(0.1)
end
pause(1)

% create expt struct
fprintf('Initializing expt struct... ')
expt.date = daqStartDate;
expt.time = daqStartTime;
expt.setup = setup;
expt.inputs = inputs;
expt.outputs = outputs;
fprintf('OK.\n')