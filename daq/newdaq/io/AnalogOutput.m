classdef AnalogOutput < signal
    properties
        type = 'analog';
        dir = 'output';
    end
    methods
        function obj = AnalogOutput(name, line)
            obj.name = name;
            obj.line = line;
            obj.channel = ['ao' num2str(line)];
        end
        
        function addToDaq(obj, dq)
            ch = addoutput(dq, obj.dev, obj.line, 'Voltage');
            ch.TerminalConfig = 'SingleEnded';
        end
    end
end