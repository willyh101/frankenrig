% grating settings
ori = 30;
sz = 30;
contrast = 100;
incl_blank = 1;

% x_loc = 0;
% y_loc = 0;

tf = 1;
sf = 0.08;

% other settings
dist_screen_cm = 10;
vert_screen_sz = 15;
gf = 5;
gtype = 'box';
method = 'symmetric';

% save locations
savefolder = '';
remotesavefolder = '';

% etc file paths
gammafile = '/Users/miscFrankenRig/Documents/HBIO/gamma_correction_210330.mat';

d = daq('ni');

% outputs
addoutput(d, 'Dev2', 'Port0/Line1', 'OutputOnly'); % stim on indicator to master
addoutput(d, 'Dev2', 'Port0/Line4', 'OutputOnly'); % stim on indicator to SI
% addoutput(d, 'Dev2', 'Port0/Line5', 'OutputOnly'); % unused trigger line to master
addoutput(d, 'Dev2', 'Port0/Line3', 'OutputOnly'); % pulsing info to master

% inputs
addinput(d, 'Dev2', 2, 'Voltage'); % unsure of fxn
addtrigger(d, 'Digital', 'StartTriggers', 'External', 'Dev2/PFI1'); % trigger channel
% d.NumDigitalTriggersPerRun = 100; % why 100?

%% get wininfo and etc
AssertOpenGL;
screens = Screen('Screens');
screenNumber = max(screens);

Screen('Preference', 'SkipSyncTests', 1);

% Define black and white
white = WhiteIndex(screenNumber);
black = BlackIndex(screenNumber);
gray=round((white+black)/2);

if gray == white
    gray=white / 2;
end

inc=white-gray;

% open a window
[w,w_rect] = PsychImaging('OpenWindow',screenNumber,gray);


% get window size
xRes = RectWidth(Screen('Rect', screenNumber));
yRes = RectHeight(Screen('Rect', screenNumber));
vertScreenDimDeg = atand(vert_screen_sz/dist_screen_cm);
pixPerDeg = yRes/vertScreenDimDeg;

fr = Screen('FrameRate', screenNumber);
ifi = Screen('GetFlipInterval', w);
Bcol = 128;
bg = ones(yRes,xRes)*Bcol;
bg = Screen('MakeTexture', w, bg);

% load gamma
% load(gammaFile,'gamma_table')
% Screen('LoadNormalizedGammaTable',w,gamma_table);


%% make grating

tex = make_grating(wininfo, grating)








%%make combos (square)

patchRadiusPx = ceil(sz*pixPerDeg/2);
numFrames = ceil(fr/tf);

% assume it's the center of the screen for now
[x0, y0] = RectCenter(w_rect);

clear tex
for i=1:numFrames
    phase = (i/numFrames)*2*pi;
    [x,y] = meshgrid(-patchRadiusPx:patchRadiusPx,-patchRadiusPx:patchRadiusPx);
    angle=ori*pi/180;
    f=sf*2*pi; % cycles/pixel
    a=cos(angle)*f;
    b=sin(angle)*f;
%     m = exp(-((x/(gf*patchRadiusPx)).^2)-((y/(gf*patchRadiusPx)).^2));
    m=exp(-((x/90).^2)-((y/90).^2)).*sin(a*x+b*y+phase);
    tex(i)=Screen('MakeTexture', w, round(gray+inc*m)); 
%     if streq(gtype, 'box')
%         s = sin(a*x+b*y+phase);
%         ext = max(max(max(s)),abs(min(min(s))));
%         g0=ext*((s>0)-(s<0));%.*g0;
%     end
%     if streq(method, 'symmetric')
%         incmax = min(255-Bcol,Bcol);
%         g0 = (floor(contrast*(incmax*g0)+Bcol));
%     end
    tex(i) = Screen('MakeTexture', w, round(gray+inc*m));
end


%%alternate way
% p=ceil(1/tf);
% frad=tf*2*pi;
% visiblesize=2*patchRadiusPx+1;
% x = meshgrid(-patchRadiusPx:patchRadiusPx + p, 1);
% grating=gray + 255*cos(frad*x);
% 
% 
% gratingtex=Screen('MakeTexture', w, grating);

movieDurationSecs=5;
frameRate=Screen('FrameRate',screenNumber);

movieDurationFrames=round(movieDurationSecs * frameRate);
movieFrameIndices=mod(0:(movieDurationFrames-1), numFrames) + 1;

% Use realtime priority for better timing precision:
    priorityLevel=MaxPriority(w);
    Priority(priorityLevel);

    % Animation loop:
    for i=1:movieDurationFrames
        % Draw image:
        Screen('DrawTexture', w, tex(movieFrameIndices(i)));
        % Show it at next display vertical retrace. Please check DriftDemo2
        % and later, as well as DriftWaitDemo for much better approaches to
        % guarantee a robust and constant animation display timing! This is
        % very basic and not best practice!
        Screen('Flip', w);
    end

    Priority(0);

    % Close all textures. This is not strictly needed, as
    % Screen('CloseAll') would do it anyway. However, it avoids warnings by
    % Psychtoolbox about unclosed textures. The warnings trigger if more
    % than 10 textures are open at invocation of Screen('CloseAll') and we
    % have 12 textues here:
    Screen('Close');

    % Close window:
    sca;



















