t = 1;
hz = 40;
minsep = 7;

spks = poissonSpikeGenerator(t,hz,minsep);

disp(['Mean rate (Hz): ' num2str(length(spks)/t)])
