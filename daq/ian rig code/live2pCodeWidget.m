load('T:\live2p_out\data.mat')

out.PSTHs = onlinePSTHs;
out.onlineTraces = onlineTraces;
out.onlineTraceLengths = onlineTrialLengths
%%
ExpStruct.orientationData50 = out;

%%
ExpStruct.spikeTest = out

%%
ExpStruct.ManifoldWrite2 = out