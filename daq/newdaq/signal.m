classdef signal < handle
    
    properties (Abstract)
        type
        dir
    end

    properties
        name
        dev
        port
        line
        channel
        notes
    end

    methods
        function obj = signal(name, dev, port, line)
            obj.name = name;
            obj.dev = dev;
            obj.port = port;
            obj.line = line;
            obj.channel = ['port' num2str(port) '/line' num2str(line)];
        end
    end

end