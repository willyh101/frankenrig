classdef AnalogPulseOuput < PulseOutput

    properties (Dependent)
        pulse_lengths_samples
    end

    methods
        function sam = get.pulse_lengths_samples(obj)
            sam = obj.to_samples(obj.pulse_starts);
        end

        function addConstantPulse(obj, start, length, ampl)
            pulse_start_samples = obj.to_samples(start);
            pulse_length_samples = obj.to_samples(length);
            obj.pulse_starts = [obj.pulse_starts, start];
            obj.sweep(pulse_start_samples:pulse_start_samples+pulse_length_samples-1) = ampl;
        end

        function addVaryingPulse(obj, start, vals)
            pulse_start_samples = obj.to_samples(start);
            length = numel(vals);
            pulse_length_samples = obj.to_samples(length);
            obj.pulse_starts = [obj.pulse_starts, start];
            obj.sweep(pulse_start_samples:pulse_start_samples+pulse_length_samples-1) = vals;
        end
    end
end