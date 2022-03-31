function powers = autoCalibSIPowerFun(startPower,zPowerReference,z,lz)
global autoCalibPlaneToUse
if isempty(autoCalibPlaneToUse)
    autoCalibPlaneToUse =0;
end

% disp(['autoCalibPower Ran z: ' num2str(autoCalibPlaneToUse)]);

powers = ones(size(z))*0.1;

try
powers(z==autoCalibPlaneToUse)=startPower;
catch
end

% disp(num2str(powers))
% powers = sqrt(z-zPowerReference).*startPower;



end