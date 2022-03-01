function streamFrameData(src,evt,varargin)

hSI = src.hSI;
lastStripe = hSI.hDisplay.stripeDataBuffer{hSI.hDisplay.stripeDataBufferPointer};
frame = lastStripe.roiData{1}.imageData{1}{1};
z = lastStripe.roiData{1}.zs;
f = lastStripe.roiData{1}.frameNumberAcq;
% figure(1111)
% clf
% imagesc(frame)
% title(z)
% ylabel(f)
% drawnow