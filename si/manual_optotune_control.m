hSI.abort();
hSI.hBeams.pzAdjust = 1;
hSI.hBeams.pzCustom= {@autoCalibSIPowerFun} ;

global autoCalibPlaneToUse

autoCalibPlaneToUse = 60;

hSI.hStackManager.arbitraryZs = [ 0 autoCalibPlaneToUse];

hSI.startGrab();