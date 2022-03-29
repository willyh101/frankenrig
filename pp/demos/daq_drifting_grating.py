from psychopy import visual, core, event
from psychopy.visual import grating
from psychopy.monitors import Monitor
import nidaqmx as ni
from nidaqmx.constants import Signal


stim_time = 1
tf = 1
pre_time = 1

mon = Monitor('iPad Retina')
mon.setSizePix([1600,1200])
mon.setDistance(10)
mon.setWidth(19.7)

win = visual.Window(size=[1920,1080], units='deg', screen=0, allowGUI=False, monitor=mon)

stim = grating.GratingStim(
    win, 
    tex='sqr',
    contrast=1,
    ori=0,
    size=30,
    sf=0.08
)

clock = core.Clock()

def display_grating():
    win.flip()
    
    clock.reset()
    while clock.getTime() < pre_time:
        pass
    
    clock.reset()
    while clock.getTime() < stim_time:
        t = clock.getTime()
        stim.phase = tf*t
        stim.draw()
        win.flip()
        
    win.flip()

# while True:
#     for keys in event.getKeys():
#         if keys in ['escape', 'q']:
#             core.quit()
            
with ni.Task() as task:
    task = ni.Task()
    task.di_channels.add_di_chan('Dev....')
    task.register_signal_event(Signal.START_TRIGGER, display_grating)
    task.read()

win.close()