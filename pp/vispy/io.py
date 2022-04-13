# Digital Outputs
SI_FLIP = 'Dev2/Port0/Line4'
DAQ_FLIP = 'Dev2/Port0/Line1'
SI_TRIG = 'Dev2/Port0/Line2'
DAQ_PULSE = 'Dev2/Port0/Line3'
DAQ_CLOCK = 'Dev2/Port0/Line0'

# Triggers
PT_TRIG = {
    'pfi': '/Dev2/PFI1',
    'ai': 'Dev2/ai0'
}

def make_standard_config():
    digital_outs = [SI_FLIP, DAQ_FLIP, DAQ_CLOCK, DAQ_PULSE]
    flip = [1, 1, 0, 0]
    blank = [0, 0, 0, 0]
    pulse_on = [0, 0, 1, 1]
    pulse_off = [0, 0, 1, 0]
    return (digital_outs, flip, blank, pulse_on, pulse_off)