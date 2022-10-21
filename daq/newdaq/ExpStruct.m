classdef ExpStruct < handle
    properties
        date
        time
        setup
        daq
        inputs
        outputs
        data
        notes
    end
    
    methods
        function obj = ExpStruct(setup, inputs, outputs)
            fprintf('Initializing expt struct... ')
            exp_datetime = datetime('now');
            obj.date = datestr(exp_datetime, 'mm/dd/yyyy');
            obj.time = datestr(exp_datetime, 'HH:MM:SS');
            obj.setup = setup;
            obj.inputs = inputs;
            obj.outputs = outputs;
            fprintf('OK.\n')
        end
       
    end
end