from psychopy.visual import grating
import nidaqmx as ni
import numpy as np

import sys
sys.path.append('m:/will/rigcode/pp')

from vispy.monitors import ipad_monitor
from vispy.utils import display_blank_trial, display_grating, make_all_combos

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

# trigger pt (in task)
# NI-6003
ai = 'Dev2/ai0'
pfi = '/Dev2/PFI1'
sample_rate = 1000
timeout = 30

# notify daq (out task)
digital_outs = {
    'si_flip': 'Dev2/Port0/Line4',
    'daq_flip': 'Dev2/Port0/Line1',
    'daq_trig': 'Dev2/Port0/Line5',
    'si_trig': 'Dev2/Port0/Line2'
}

# make generic output signals
flip_on = [True, True, False, False]
flip_off = [False, False, False, False]

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
    trig_task.ai_channels.add_ai_voltage_chan(ai)
    trig_task.triggers.start_trigger.cfg_dig_edge_start_trig(pfi)
    trig_task.timing.cfg_samp_clk_timing(sample_rate, samps_per_chan=2)

    # setup digital outputs
    for do in digital_outs.values():
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
            out_task.write(flip_off)
        
        # flip a grey screen
        win.flip()

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

# then save the stim record into frankenlocal/tmp