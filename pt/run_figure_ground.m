clear

mouse = 'w46_4';
date = '20220518';
notes = 'figure ground with hole, 16 conditions, iso and cpx, fg3, hopefully still aligned ret and imaging';

fsize = 15;
gsize = 50;

oris = 0:45:135;

sf = 0.08;
tf = 1;

vertScreenSize = 15;
dScreen = 10;

duration = 1;
pre_wait = 1;
post_wait = 3;

maxNumTrials = 650;

screen = 1;

result.mouse = mouse;
result.date = date;
result.notes = notes;
result.fsize = fsize;
result.gsize = fsize;
result.oris = oris;
result.sf = sf;
result.tf = tf;
result.vertScreenSize = vertScreenSize;
result.dScreen = dScreen;
result.duration = duration;
result.pre_wait = pre_wait;
result.post_wait = post_wait;


% daq stuff
dq = daq.createSession('ni');
addDigitalChannel(dq,'Dev2', 'Port0/Line1', 'OutputOnly'); %stim on indicator to master
addDigitalChannel(dq,'Dev2', 'Port0/Line4', 'OutputOnly'); %stim on indicator to SI
addDigitalChannel(dq,'Dev2', 'Port0/Line5', 'OutputOnly'); %unused trig line to master
addDigitalChannel(dq,'Dev2', 'Port0/Line3', 'OutputOnly'); % pulsing info to master
addDigitalChannel(dq,'Dev2', 'Port0/Line2', 'OutputOnly'); %trig to SI

% window setup
black = BlackIndex(screen);
white = WhiteIndex(screen);
gray=round((white+black)/2);
inc=white-gray;
PsychImaging('PrepareConfiguration');
Screen('Preference', 'SkipSyncTests',1);
Screen('Preference','VisualDebugLevel',3);
xRes = RectWidth(Screen('Rect', screen));
yRes = RectHeight(Screen('Rect', screen));
[w, screenRect]=Screen('OpenWindow',screen, gray);
vScreenDeg = atand(vertScreenSize/dScreen);
pxPerDeg = yRes/vScreenDeg;
frameRate = Screen('FrameRate',screen);
gammaFile = 'C:/Users/miscFrankenRig/Documents/HBIO/gamma_correction_210330.mat';
load(gammaFile,'gamma_table')
Screen('LoadNormalizedGammaTable',w,gamma_table);
AssertOpenGL;
Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

x0 = floor(xRes/2);
y0 = floor(xRes/2);

mask_radius = ceil(fsize*pxPerDeg/2);
outer_mask_radius = ceil(gsize*pxPerDeg/2);

bg = ones(yRes,xRes)*gray;
bg = Screen('MakeTexture', w, bg);

% Calculate parameters of the Sinusoid:
sfPxPerCyc = 1/sf*pxPerDeg;
p = ceil(sfPxPerCyc);
fr = 1/p*2*pi; 
visiblesize=2*mask_radius+1;
outer_visiblesize = 2*outer_mask_radius+1;

% Create one single static Sinusoid image:
[x,y]=meshgrid(-mask_radius:mask_radius + p, -mask_radius:mask_radius);
sinusoid=gray + inc*cos(fr*x);

% Create surround grating
[x,y]=meshgrid(-outer_mask_radius:outer_mask_radius + p, -outer_mask_radius:outer_mask_radius);
outer_sinusoid=gray + inc*cos(fr*x);

grating_mask = sinusoid>gray; %returns a 0 and 1 mask
outer_grating_mask = outer_sinusoid>gray; %returns a 0 and 1 mask

contrast_diff = 1;
sinusoid(grating_mask) = gray-(255*contrast_diff/2); %makes actual squarewave
sinusoid(~grating_mask) = gray+(255*contrast_diff/2);
outer_sinusoid(outer_grating_mask) = gray-(255*contrast_diff/2); %makes actual squarewave
outer_sinusoid(~outer_grating_mask) = gray+(255*contrast_diff/2);


sinusoidtex=Screen('MakeTexture', w, sinusoid);
% make gray center for inverse aperture
gray_center = ones(size(sinusoid));
gray_center = gray_center*gray;
gray_center_texture = Screen('MakeTexture', w, gray_center);

%make circular transparency mask texture
%     masktex=Screen('MakeTexture', w, mask);
masktex=Screen('MakeTexture', w, outer_sinusoid);

x_shift = 0;
y_shift = 0;
% Definition of the drawn rectangle on the screen:
dstRect=[0 0 visiblesize visiblesize];
dstRect=CenterRect(dstRect, screenRect);
[x_center y_center] = RectCenterd(dstRect);
dstRect = CenterRectOnPointd(dstRect,x_center+x_shift,y_center+y_shift); %%%%%%%

% make centered gray patch for inverse aperture 20% larger than CRF
% mask
dstRectGrey=dstRect*1.5;
dstRectGrey=CenterRect(dstRectGrey, screenRect);
[x_centerG y_centerG] = RectCenterd(dstRectGrey);
dstRectGrey = CenterRectOnPointd(dstRectGrey,x_centerG+x_shift,y_centerG+y_shift); %%%%%%%


outer_dstRect=[0 0 outer_visiblesize outer_visiblesize];
outer_dstRect=CenterRect(outer_dstRect, screenRect);
[x_center y_center] = RectCenterd(outer_dstRect);
outer_dstRect = CenterRectOnPointd(outer_dstRect,x_center+x_shift,y_center+y_shift); %%%%%%%

% Query duration of monitor refresh interval:
ifi=Screen('GetFlipInterval', w);

shiftperframe= tf * p * ifi;

% wait to make sure everything is synched
WaitSecs(5);

for tr=1:maxNumTrials
    stimID = randi(4);
    
    if numel(oris) > 1
        direction = randsample(oris,1);
    else
        direction = oris;
    end
    
    result.stims(tr) = stimID;
    result.dirs(tr) = direction;

    disp(['Trial ' num2str(tr)])
    disp(['Stim ID ' num2str(stimID)])
    
    outputSingleScan(dq,[0 0 0 0 1])
    outputSingleScan(dq,[0 0 0 0 0])

    
    i = 0;

    % Perform initial Flip to sync us to the VBL and for getting an initial
    % VBL-Timestamp for our "WaitBlanking" emulation:
    
    WaitSecs(pre_wait);

    % We run at most 'movieDurationSecs' seconds if user doesn't abort via keypress.
    outputSingleScan(dq,[1 1 0 0 0])
    vbl=Screen('Flip', w);
    vblendtime = vbl + duration;
    
    while (vbl<vblendtime)
        xoffset = mod(i*shiftperframe,p);
        i=i+1;

        srcRect=[xoffset 0 xoffset + visiblesize visiblesize];
        outer_srcRect=[xoffset 0 xoffset + outer_visiblesize outer_visiblesize];
        switch stimID
            case 1
                % iso mask
                Screen('DrawTexture', w, masktex, outer_srcRect, outer_dstRect, direction);

            case 2
                % cross mask
                Screen('DrawTexture', w, masktex, outer_srcRect, outer_dstRect, direction+90);
                Screen('DrawTexture', w, sinusoidtex, srcRect, dstRect, direction);
                

            case 3
                % center
                Screen('DrawTexture', w, sinusoidtex, srcRect, dstRect, direction);

            case 4
                % hole
                Screen('DrawTexture', w, masktex, outer_srcRect, outer_dstRect, direction);
                Screen('DrawTexture', w, gray_center_texture, srcRect, dstRect, direction);
        end
        
        vbl=Screen('Flip', w);
        

    end
    Screen('FillRect', w, gray, []);
    Screen('Flip',w);
    outputSingleScan(dq,[0 0 0 0 0])

    WaitSecs(post_wait);

end

sca