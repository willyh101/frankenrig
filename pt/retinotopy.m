% grating settings
ori = 30; % [0 30 60];
sz = 30; %; [30, 50];
tf = 1;
sf = 0.08;
contrast = 100; % [100];

% save settings

% setup daq comm
d = daq('ni');

% output
addoutput(d, 'Dev2', 'Port0/Line1', 'OutputOnly'); % stim on indicator to master
addoutput(d, 'Dev2', 'Port0/Line4', 'OutputOnly'); % stim on indicator to SI
addoutput(d, 'Dev2', 'Port0/Line5', 'OutputOnly'); % unused trigger line to master
addoutput(d, 'Dev2', 'Port0/Line3', 'OutputOnly'); % pulsing info to master

% input
addinput(d, 'Dev2', 2, 'Voltage'); % unsure of fxn
addtrigger(d, 'Digital', 'StartTriggers', 'External', 'Dev2/PFI1'); % trigger channel

d.NumDigitalTriggersPerRun = 100; % why 100?
% s0.TriggersPerRun = 100;
% s0.ExternalTriggerTimeout = 120;
% s0.NumberOfScans=2;


% load gamma

%% setup stims
conds = {ori, sz, tf, sf, contrast};
nvars = numel(conds);
neach = cellfun(@numel, conds);
ncombo = prod(neach);
combos = zeros(nvars,ncombo);

for j=1:ncombo
    for i=1:nvars
        ind = mod(floor((j-1)/prod([1 neach(1:i-1)])),neach(i))+1;
        combos(i,j) = conds{i}(ind);
    end
end

[x,y] = meshgrid(-300:300, -300:300);
nf = 100;
x=1;
y=1;
gf=0;
bg = 172*ones(1000,1000);

clear T G

for i=1:nf
    phase = (i/nf)*2*pi;
    angle = ori*pi/180;
    f = sf/10*2*pi;
    a = cos(angle)*f;
    b = sin(angle)*f;
    g0 = exp(-((x/(gf*300)).^2)-((y/(gf*300)).^2));
    s = sin(a*x+b*y+phase);
    ext = max(max(max(s)),abs(min(min(s))));
    G0=ext*((s>0)-(s<0));
    incmax = min(255-172,172);
    G = (floor(contrast*(incmax*G0)+172));
    T = bg;
    T(1:1+size(G,2)-1, 1:1+size(G,2)-1) = G;
end



% PT try-catch
