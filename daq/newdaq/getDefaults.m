function s = getDefaults()

% nidaq configs
s.daqrate = 20000;
s.device = 'Dev1';

% paths
s.power_calib_path = 'path/to/power_calibration';
s.spatial_calib_path = 'path/to/spatial_calibration';
s.holorequest_path = 'path/to/holorequest';
s.save_path = 'path/to/save/results';

% record keeping
s.rig = 'frankenscope';
s.stim_laser = 'monaco 40 W';

% laser and other constants
s.eom_offset = -0.15;