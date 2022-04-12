import nidaqmx as ni
import time

digital_outs = {
    'si_flip': 'Dev2/Port0/Line4',
    'daq_flip': 'Dev2/Port0/Line1',
    # 'daq_trig': 'Dev2/Port0/Line5',
    # 'si_trig': 'Dev2/Port0/Line2',
    'daq_pulse': 'Dev2/Port0/Line3'
}

# make generic output signals
clear = [False, False, False]
flip_on = [True, True, False]
flip_off = [False, False, False]
pulse_on = [False, False, True]
pulse_off = [[False, False, False, False, False, False],
             [False, False, True, False, True, False],
             [False, False, False, False, False, False]]


with ni.Task() as task:
    for do in digital_outs.values():
        task.do_channels.add_do_chan(do)
    task.start()

    print(task.do_channels.channel_names)
    input('starting test. press enter when ready.')
    task.write(clear)

    time.sleep(2)
    print('sending flip on')
    task.write(flip_on)
    
    time.sleep(2)
    print('sending flip off')
    task.write(flip_off)

    time.sleep(2)
    print('sending 2 daq pulses')
    task.write(pulse_off)

    time.sleep(2)
    task.write(clear)

print('done.')