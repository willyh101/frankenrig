classdef DAQoutput < DAQgroup

    methods
        function obj = DAQoutput(dq, name)
            obj = obj@DAQgroup(dq)
            if nargin > 1
                obj.group = name;
            end
        end

        function tg = createPulseOutput(obj, len)
            tg = PulseOutput(obj.dq.Rate, len);
        end

        function swp = createSweep(obj, len)
            swp = zeros(len*obj.dq.Rate, obj.nch);
        end

    end
end