from psychopy import visual, core, event
from psychopy.visual import grating
from psychopy.monitors import Monitor

stim_time = 1
tf = 1
pre_time = 1

mon = Monitor('testMon')
mon.setSizePix([1920,1080])
mon.setDistance(10)
mon.setWidth(15)

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
clock.reset()

while clock.getTime() < pre_time:
    win.flip()
    
clock.reset()
while clock.getTime() < stim_time:

    t = clock.getTime()
    stim.phase = tf*t
    stim.draw()
    win.flip()
    
    for keys in event.getKeys():
        if keys in ['escape', 'q']:
            core.quit()
            
win.close()