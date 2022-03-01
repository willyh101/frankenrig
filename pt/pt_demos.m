%% Demo 1
% basic guide on how to display animations of a drifting grating

try
    AssertOpenGL;
    screens=Screen('Screens'); % screen 0 is the one with the menu bar
    screenNumber=max(screens);

    % find white and black vals
    white = WhiteIndex(screenNumber);
    black = BlackIndex(screenNumber);
    gray=round((white+black)/2);
    
    % This makes sure that on floating point framebuffers we still get a
    % well defined gray. It isn't strictly neccessary in this demo:
    if gray == white
        gray=white / 2;
    end
    
    % Contrast 'inc'rement range for given white and gray values:
    inc=white-gray;
    
    % opem a double buffered window with grey background
%     w=Screen('OpenWindow',screenNumber, gray);

    % Compute each frame of the movie and convert the those frames, stored in
    % MATLAB matices, into Psychtoolbox OpenGL textures using 'MakeTexture';
    figure(1);
    clf

    numFrames=12; % temporal period, in frames, of the drifting grating
    for i=1:numFrames
        clf
        phase=(i/numFrames)*2*pi;
        % grating
        [x,y]=meshgrid(-300:300,-300:300);
        angle=30*pi/180; % 30 deg orientation.
        f=0.05*2*pi; % cycles/pixel
        a=cos(angle)*f;
        b=sin(angle)*f;
        m=exp(-((x/90).^2)-((y/90).^2)).*sin(a*x+b*y+phase);
        tex = round(gray+inc*m);
        image(tex)
        pause

%         tex(i)=Screen('MakeTexture', w, round(gray+inc*m));
    end

catch
    sca;
    Priority(0);
    psychrethrow(psychlasterror);
end

%% Demo 2
% draw 2d grating online by using 1d grating texture

gratingsize = 400;
drawmask = 1;
f = 0.05; % cycle/pixel
cyclespersecond = 1;
angle = 30;
movieDurationSecs = 20;
texsize = gratingsize/2;

try
    AssertOpenGL;

    screens=Screen('Screens'); % screen 0 is the one with the menu bar
    screenNumber=max(screens);

    % find white and black vals
    white = WhiteIndex(screenNumber);
    black = BlackIndex(screenNumber);
    gray=round((white+black)/2);
    
    % This makes sure that on floating point framebuffers we still get a
    % well defined gray. It isn't strictly neccessary in this demo:
    if gray == white
        gray=white / 2;
    end
    
    % Contrast 'inc'rement range for given white and gray values:
    inc=white-gray;

    % First we compute pixels per cycle, rounded up to full pixels, as we
    % need this to create a grating of proper size below:
    p=ceil(1/f);

    fr=f*2*pi;
    
    % This is the visible size of the grating. It is twice the half-width
    % of the texture plus one pixel to make sure it has an odd number of
    % pixels and is therefore symmetric around the center of the texture:
    visiblesize=2*texsize+1;

    % Create one single static grating image:
    %
    % We only need a texture with a single row of pixels(i.e. 1 pixel in height) to
    % define the whole grating! If the 'srcRect' in the 'Drawtexture' call
    % below is "higher" than that (i.e. visibleSize >> 1), the GPU will
    % automatically replicate pixel rows. This 1 pixel height saves memory
    % and memory bandwith, ie. it is potentially faster on some GPUs.
    %
    % However it does need 2 * texsize + p columns, i.e. the visible size
    % of the grating extended by the length of 1 period (repetition) of the
    % sine-wave in pixels 'p':
    x = meshgrid(-texsize:texsize + p, 1);

    % Compute actual cosine grating:
    grating=gray + inc*cos(fr*x);
    
catch
    sca;
    Priority(0);
    psychrethrow(psychlasterror);
end


%% Demo 6
% surround + center stimulus


