classdef signal < handle
    
    properties (Abstract)
        type
        dir
    end
    
    properties
        name
        dev
        line
        idx
        channel
        notes
    end

    methods
        function save(obj, file_path)
            warning('off', 'MATLAB:structOnObject')
            out = struct(obj);
            warning('on', 'MATLAB:structOnObject')
            save(file_path, out, '-v7')
        end
        
        function d = get.dev(obj)
            if isempty(obj.dev)
                d = 'Dev1';
            else
                d = obj.dev;
            end
        end

        function idx = subsindex(obj)
            idx = obj.idx-1;
        end
    end
end