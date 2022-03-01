function trialLengths(src,evt,varargin)

hSI = src.hSI;
lastStripe = hSI.hDisplay.stripeDataBuffer{hSI.hDisplay.stripeDataBufferPointer};
f = lastStripe.roiData{1}.frameNumberAcq;
nplanes = numel(hSI.hStackManager.arbitraryZs);

f = f/nplanes;

savepath = hSI.hScan2D.logFilePath;
savestem = hSI.hScan2D.logFileStem;
savename = [savepath '\'  savestem '_fileLengths.csv'];

dlmwrite(savename, f, 'delimiter', ',', '-append')