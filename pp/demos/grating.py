from psychopy import visual, core
from psychopy.visual import grating

win = visual.Window(size=[1920,1080], units='pix', screen=0)
stim = grating.GratingStim(
    win, 
    tex='sqr',
    contrast=1,
    ori=35,
    size=500
)

stim.draw()
win.flip()
core.wait(5.0)