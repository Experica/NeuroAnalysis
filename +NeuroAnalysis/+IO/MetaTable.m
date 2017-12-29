classdef MetaTable < handle
    %METATABLE Metadata table
    %   Detailed explanation goes here
    
    properties (Access = private)
        Tests       % Structure of metadata
        map         % key-value map of filenames to indices
    end
    
    properties (Access = public)
       filepath     % where the metadata file should be stored
    end
    
    methods (Access = public)
        
        function obj = MetaTable(filepath)
            %METATABLE Create a MetaTable object from the given path
            % If the path does not point to a file, then an empty object
            % will be created
            obj.filepath = filepath;
            if exist(obj.filepath, 'file')
                metafile = load(obj.filepath);
                obj.Tests = metafile.Tests;
                obj.map = metafile.map;
            else
                obj.Tests = struct;
                obj.Tests.files = {};
                obj.map = containers.Map('KeyType', 'char', 'ValueType', 'double');
            end
        end
        
        function addEntry(obj, test)
            %ADDENTRY Add a test structure with any name value pair fields

            % Search for test in existing table
            fields = setdiff(fieldnames(obj.Tests), fieldnames(test));
            if obj.map.isKey(test.filepath)
                idx = obj.map(test.filepath);
                % Copy any columns that aren't present in this test
                for f=1:length(fields)
                    test.(fields{f}) = obj.Tests.(fields{f})(idx);
                end
            else
                idx = length(obj.map) + 1;
                % Fill in any columns that aren't present in this test
                for f=1:length(fields)
                    test.(fields{f}) = missing;
                end
            end
            obj.map(test.filepath) = idx;            
            
            % Change characters to strings
            fields = fieldnames(test);
            for f=1:length(fields)
                if ischar(test.(fields{f}))
                    test.(fields{f}) = string(test.(fields{f}));
                end
            end

            % Add the new test
            if length(obj.Tests.files) < idx || isempty(obj.Tests.files(idx))
                obj.Tests.files{idx} = test.filepath;
            else
                obj.Tests.files{idx} = union(obj.Tests.files{idx}, test.filepath);
            end
            test = rmfield(test, {'filepath', 'files'});
            fields = fieldnames(test);
            for f=1:length(fields)
                obj.Tests.(fields{f})(idx) = test.(fields{f});
            end

        end
        
        function Tests = query(obj, varargin)
            %QUERY Find matching tests
            nargs = length(varargin);
            if round(nargs/2)~=nargs/2
               error('Query needs propertyName/propertyValue pairs')
            end
            keys = varargin(1:2:end);
            values = varargin(2:2:end);
            isMatch = ones(1,length(obj.map));
            for kv = 1:length(keys)
                isMatch = isMatch & ismember(obj.getfieldi(obj.Tests,keys{kv}), values{kv});
            end
            Tests = structfun(@(x)x(isMatch), obj.Tests, 'UniformOutput', false);
        end
        
        function export(obj)
            %EXPORT save to disk
            disp(['Saving metadata:    ',obj.filepath,'    ...']);
            metafile.Tests = obj.Tests;
            metafile.map = obj.map;
            save(obj.filepath, '-struct', 'metafile', '-v7.3');
            disp('Saving metadata:    Done.');
        end
    end
    
    methods (Static, Access = private)
        
        function value = getfieldi(S,field)
            %GETFIELDI case-insensitive dynamic structure field access
            names = fieldnames(S);
            isField = strcmpi(field,names);  

            if any(isField)
                assert(sum(isField) == 1);
                value = S.(names{isField});
            else
                value = [];
            end
        end
        
    end
    
end

