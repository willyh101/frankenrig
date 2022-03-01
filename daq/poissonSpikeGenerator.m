function spikeTimes = poissonSpikeGenerator(hz, t, minsep, plotMe)
% generate poisson spike train with mean firing rate
% hz = mean firing rate in Hz
% t = length of stim (seconds)
% minsep = minimum seperation of spikes (ms)
% returns spike times in ms
if nargin < 4
    plotMe = 0;
end

c = 0;
while 1
    c = c+1;
    if c > 100000
        error('Could not find a spike train with a suitable min seperation!')
    end
    dt = 1/20000;
    nBins = floor(t/dt);
    spikeMat = rand(1, nBins) < hz*dt;

    tVec = 0:dt:t-dt;
    spikeTimes = tVec(spikeMat);
    spikeTimes = round(spikeTimes*1000);
    
    if ~any(diff(spikeTimes) < minsep)
        break
    end
end

spikeTimes = spikeTimes - min(spikeTimes);

if plotMe
    figure(55); clf
    hold all

    for spikeCount = 1:length(spikeTimes)
        plot([spikeTimes(spikeCount) spikeTimes(spikeCount)], [1-0.4 1+0.4])
    end
end