classdef PulseOutput < handle
    properties
        daq_rate
        sweep_length
        sweep
        sweep_length_samples
        pulse_starts
        pulse_starts_samples
    end

    properties (GetAccess = private)
        trig_length = 5; % ms
        trig_length_samples
    end

    methods
        function obj = PulseOutput(daq_rate, sweep_length)
            obj.daq_rate = daq_rate;
            obj.sweep_length = sweep_length;
            obj.sweep_length_samples = obj.daq_rate*sweep_length;
            obj.sweep = zeros(1,daq_rate*sweep_length);
        end

        function samples = to_samples(obj, ms)
            samples = ms * obj.daq_rate / 1000;
        end

        function tls = get.trig_length_samples(obj)
            tls = obj.daq_rate*obj.trig_length/1000;
        end

        function tss = get.pulse_starts_samples(obj)
            tss = obj.pulse_starts*obj.daq_rate./1000;
        end

        function resetOutput(obj)
            obj.sweep = zeros(1,obj.daq_rate*obj.sweep_length);
        end

        function show(obj, fig)
            if nargin < 2
                fig=5;
            end
            figure(fig)
            xs = 1:obj.sweep_length_samples;
            xs = xs/obj.daq_rate;
            plot(xs, obj.sweep)
        end

        function d = double(obj)
            d = obj.sweep;
        end
    end

    methods (Access = private)
        function addPulse(obj, trig_time)
            trig_start_samples = trig_time * obj.daq_rate / 1000;
            obj.pulse_starts = [obj.pulse_starts trig_time];
            obj.sweep(trig_start_samples:trig_start_samples+obj.trig_length_samples-1) = 1;
        end

        function addPulseTrain(obj, trig_times)
            for tt=trig_times
                obj.addPulse(tt)
            end
        end

        function addStartTrigger(obj)
            obj.sweep(1:obj.trig_length_samples) = 1;
            obj.pulse_starts = [obj.pulse_starts 1/obj.daq_rate];
        end
    end
end