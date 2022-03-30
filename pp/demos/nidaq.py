import nidaqmx as ni
from nidaqmx.constants import Signal

def test_callback(task_handle, signal_type, callback_data):
    print('hi there')
    return 0

task = ni.Task()
task.ai_channels.add_ai_voltage_chan('Dev2/ai0')
task.triggers.start_trigger.cfg_dig_edge_start_trig('/Dev2/PFI1')
task.timing.cfg_samp_clk_timing(rate=1000)
task.start()
task.register_signal_event(Signal.START_TRIGGER, test_callback)
task.read()
print('did you get it?')
task.close()