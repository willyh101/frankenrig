from psychopy.monitors import Monitor

def test_monitor():
    mon = Monitor('testMon')
    mon.setSizePix([1920,1080])
    mon.setDistance(10)
    mon.setWidth(15)
    return mon

def ipad_monitor():
    mon = Monitor('iPad Retina')
    mon.setSizePix([1600,1200])
    mon.setDistance(10)
    mon.setWidth(19.7)
    return mon