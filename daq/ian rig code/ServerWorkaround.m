function ServerWorkaround(ExpStruct)

invar = []; count=0;
while ~isstruct(invar) && count <120
invar = msrecv(ExpStruct.SISocket,.5);
count=count+1;
end

holoRequest = invar;

locations = FrankenScopeRigFile();
save([locations.HoloRequest 'holoRequest.mat'],'holoRequest');
save([locations.HoloRequest_DAQ 'holoRequest.mat'],'holoRequest');