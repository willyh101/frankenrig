from psychopy.core import Clock
import numpy as np
import itertools

def display_grating(stim, win, stim_time, tf):
    clock = Clock()
    clock.reset()
    while clock.getTime() < stim_time:
        t = clock.getTime()
        stim.phase = t * tf
        stim.draw()
        win.flip()
    win.flip()
    
def display_blank_trial(win, stim_time):
    clock = Clock()
    clock.reset()
    while clock.getTime() < stim_time:
        t = clock.getTime()
        win.flip()
    win.flip()
    
def make_all_combos(blanks=1, **kwargs):
    for k,v in kwargs.items():
        if not isinstance(v, (list, tuple, np.ndarray)):
            kwargs[k] = np.array([v])
            
    vals = list(itertools.product(*kwargs.values()))
    fields = list(kwargs.keys())
    conds = [dict(zip(fields, v)) for v in vals]
    
    if blanks:
        for _ in range(blanks):
            conds.append('BLANK')
        
    return conds

def int2bits(num, nbits=12):
    num_ = np.array(num, dtype=np.uint8)
    return np.unpackbits(num_, count=nbits)