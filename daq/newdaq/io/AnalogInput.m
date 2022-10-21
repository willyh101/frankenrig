classdef AnalogInput < signal
    properties
        type = 'analog';
        dir = 'input';
    end
    methods
        function obj = AnalogInput(name, line)
            obj.name = name;
            obj.line = line;
            obj.channel = ['ai' num2str(line)];
        end
        
        function addToDaq(obj, dq)
            ch = addinput(dq, obj.dev, obj.line, 'Voltage');
            ch.TerminalConfig = 'SingleEnded';
        end
    end
end