classdef EdgeCounter < signal
    properties
        type = 'counter';
        dir = 'input';
        edge = 'RisingEdge';
        ctr
        loc
    end

    methods
        function obj = EdgeCounter(name, ctr)
            obj.name = name;
            obj.ctr = ctr;
            obj.channel = 'n/a';
        end

        function addToDaq(obj, dq)
            addinput(dq, obj.dev, obj.ctr, 'EdgeCount')
        end
    end
end

