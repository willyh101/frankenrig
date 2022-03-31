import nidaqmx as ni
from nidaqmx.constants import Signal

def test_callback(task_handle, signal_type, callback_data):
    print('hi there')
    return 0

task = ni.Task()
task.ai_channels.add_ai_voltage_chan('Dev2/ai0')
task.triggers.start_trigger.cfg_dig_edge_start_trig('/Dev2/PFI1')
task.triggers.start_trigger.max_num_trigs_to_detect = 1000
task.timing.cfg_samp_clk_timing(rate=1000)
task.timing.register_done_event(test_callback)
task.start()
task.read(number_of_samples_per_channel=1)
print('did you get it?')
task.close()


# or does it need to be in a loop??
task = ni.Task()
task.ai_channels.add_ai_voltage_chan('Dev2/ai0')
task.triggers.start_trigger.cfg_dig_edge_start_trig('/Dev2/PFI1')
task.timing.register_every_n_samples_acquired_into_buffer_event(100, test_callback)
task.start()
for i in range(1000):
    task.read(number_of_samples_per_channel=100)
    
class TriggerCallback:
    def __init__(self, input_chan, trig_chan) -> None:
        self.input_chan = input_chan
        self.trig_chan = trig_chan
        self.trigger_after = 100
        self.trigger_func = None
        self.max_num_triggers = 1000
    
    def set_trigger_func(self, func):
        self.trigger_func = func
        
    def _prepare(self):
        task = ni.Task()
        task.ai_channels.add_ai_voltage_chan(self.input_chan)
        task.triggers.start_trigger.cfg_dig_edge_start_trig('/Dev2/PFI1')
        task.timing.register_every_n_samples_acquired_into_buffer_event(self.trigger_after, self.trigger_func)
        return task
    
    def run(self):
        task = self._prepare()
        task.start()
        for _ in range(self.max_num_triggers):
            task.read(number_of_samples_per_channel=self.trigger_after)
        # ideally this should be context managed
        task.stop()
        task.close()