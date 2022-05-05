classdef Trigger < signal
    properties
        type = 'digital';
        dir = 'output';
        trigTime
        trigTimeBase = 'ms';
        trigLength = 1;
    end
    
    methods
    end
end
