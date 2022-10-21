classdef Setup < handle
    properties
        daqrate = 1000;
        device = 'Dev1';
        power_calib_path = 'path/to/power_calibration';
        spatial_calib_path = 'path/to/spatial_calibration';
        holorequest_path = 'path/to/holorequest';
        save_path = 'path/to/save/results';
        rig
        stim_laser
        eom_offset
        etc
    end
end

