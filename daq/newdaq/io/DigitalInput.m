classdef DigitalInput < signal

    properties
        type = 'digital';
        dir = 'input';
        port
    end

    methods
        function obj = DigitalInput(name, port, line)
            obj.name = name;
            obj.port = port;
            obj.line = line;
            obj.channel = ['port' num2str(obj.port) '/line' num2str(obj.line)];
        end

        function addToDaq(obj, dq)
            addinput(dq, obj.dev, obj.channel, 'Digital');
        end
    end

end