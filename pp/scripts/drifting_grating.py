from psychopy.visual import grating
import nidaqmx as ni
import numpy as np

import sys
sys.path.append('m:/will/rigcode/pp')

from vispy.monitors import ipad_monitor
from vispy.utils import display_blank_trial, display_grating, make_all_combos, int2bits
from vispy.io import make_standard_config, PT_TRIG

# GENERAL SETTINGS
stim_time = 1
n_trials = 5

# STIM SETTINGS
# note: variables must match attributes of psychopy.visual.grating
ori = np.arange(0,330,30)
contrast = 1
sf = 0.08
size = 50
tf = 1
blanks = 1

# daq settings
sample_rate = 1000
timeout = 30
digital_outs, flip_on, blank, pulse_on, pulse_off = make_standard_config()

mon, win = ipad_monitor()

combos = make_all_combos(blanks=blanks, ori=ori, contrast=contrast, sf=sf, size=size)
ncombos = len(combos)
print(f'You have a total of {ncombos} conds')

# make first stim
rand_grating = np.random.choice(combos)
stim = grating.GratingStim(win, tex='sqr', **rand_grating)

# init stim record
stim_record = []

# do the main loop
# run in a context manager to avoid bad nidaq breakdowns
with ni.Task() as trig_task, ni.Task() as out_task:
    
    # setup triggers for psychopy
    trig_task.ai_channels.add_ai_voltage_chan(PT_TRIG['ai'])
    trig_task.triggers.start_trigger.cfg_dig_edge_start_trig(PT_TRIG['pfi'])
    trig_task.timing.cfg_samp_clk_timing(sample_rate, samps_per_chan=2)

    # setup digital outputs
    for do in digital_outs:
        out_task.do_channels.add_do_chan(do)
    out_task.start()

    # main loop
    tr_num=0
    for tr_num in range(n_trials):
        tr_num += 1
        print(f'trial # {tr_num} started!')
        print(f'going to display: {rand_grating}')
        win.flip()
        
        # start the nidaq
        trig_task.start()
        trig_task.read(timeout=timeout)
        
        # display grating
        if is_blank:
            display_blank_trial()
        else:
            out_task.write(flip_on)
            display_grating(stim, win)
            out_task.write(blank)
        
        # flip a grey screen
        win.flip()
        
        # stim the daq with the appropriate info
        # binvec = int2bits()

        # record the stim into the record
        stim_record.append(rand_grating)
        
        # stop the nidaq 
        trig_task.stop()
            
        # get a random grating for the next trial
        rand_grating = np.random.choice(combos)
        if isinstance(rand_grating, dict):
            for k,v in rand_grating.items():
                setattr(stim, k, v)
            is_blank = False
        else:
            is_blank = True
    
    # flip to grey screen if all the trials are done 
    win.flip()
    out_task.stop()

# then save the stim record into frankenlocal/tmp