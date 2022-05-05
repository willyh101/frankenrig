classdef DigitalInput < signal
    properties
        type = 'digital';
        dir = 'input';
    end
    methods
        % maybe move this to a exp/io group so you can add them in order
        function addToDaq(obj, dq)
            addinput(dq, obj.dev, obj.channel, 'Digital');
        end
    end
end