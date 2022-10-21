classdef DAQoutput < DAQgroup
    properties
        triggers
    end
    methods
        function obj = DAQoutput(dq, name)
            obj = obj@DAQgroup(dq)
            if nargin > 1
                obj.group = name;
            end
        end

        function tg = createPulseOutput(obj, len)
            tg = PulseOutput(obj.dq.Rate, len);
            obj.triggers{end+1} = tg;
        end

        function swp = createSweep(obj, len)
            swp = zeros(len*obj.dq.Rate, obj.nch);
        end

        function show(obj, fg)
            if nargin > 1
                figure(fg)
            else
                fg = figure(12);
            end
            hold on
            for i=1:numel(obj.triggers)
                tg = obj.triggers{i};
                tg.show(fg)
            end
        end
    end
end