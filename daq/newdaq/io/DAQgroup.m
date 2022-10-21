classdef DAQgroup < dynamicprops

    properties
        dq
        group
        names = {};
        channels = {};
    end

    properties (Dependent)
        nch
    end

    methods
        function obj = DAQgroup(dq, name)
            obj.dq = dq;
            if nargin > 1
                obj.group = name;
            end
        end

        function addio(obj, varargin)
            for i=1:numel(varargin)
                io = varargin{i};
                name = strrep(io.name, ' ', '_');
                addprop(obj, name);
                obj.(name) = io;
                obj.channels{end+1} = io;
                obj.names{end+1} = io.name;
                newIdx = numel(obj.names);
                io.idx = newIdx;
            end
        end
        
        function io_objs = getIOobjs(obj)
            io_objs = cellfun(@(x) obj.(strrep(x,' ','_')), obj.names, 'UniformOutput',false);
        end

        function addToDaq(obj)
            lines = obj.getIOobjs();
            for i=1:numel(lines)
                io = lines{i};
                io.addToDaq(obj.dq)
            end
        end    

        function idx = getIdx(obj, name)
            idx = strfind(obj.names, name);
            idx = idx{:};
        end

        function nch = get.nch(obj)
            nch = numel(obj.channels);
        end
    end
end