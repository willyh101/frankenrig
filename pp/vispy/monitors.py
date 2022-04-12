from psychopy.monitors import Monitor
from psychopy.visual import Window

def test_monitor():
    mon = Monitor('testMon')
    mon.setSizePix([1920,1080])
    mon.setDistance(10)
    mon.setWidth(15)
    return mon

def ipad_monitor(**kwargs):
    res = [1600,1200]    
    mon = Monitor('iPad Retina')
    mon.setSizePix(res)
    mon.setDistance(10)
    mon.setWidth(19.7)
    win = Window(size=res, units='deg', screen=0, allowGUI=False, monitor=mon, fullscr=True, **kwargs)
    return mon, win