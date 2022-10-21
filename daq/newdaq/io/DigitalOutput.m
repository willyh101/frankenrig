classdef DigitalOutput < signal

    properties
        type = 'digital';
        dir = 'output';
        port
        output
    end

    methods
        function obj = DigitalOutput(name, port, line)
            obj.name = name;
            obj.port = port;
            obj.line = line;
            obj.channel = ['port' num2str(obj.port) '/line' num2str(obj.line)];
        end

        function addToDaq(obj, dq)
            addoutput(dq, obj.dev, obj.channel, 'Digital');
        end
    end

end