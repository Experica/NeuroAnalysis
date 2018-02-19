classdef MetaTable < handle
    %METATABLE Metadata table
    %   Detailed explanation goes here
    
    properties (Access = private)
        Tests       % Structure of metadata
    end
    
    methods (Access = public)
        
        function obj = MetaTable(filepath)
            %METATABLE Create a MetaTable object from the given path
            % If the path does not point to a file, then an empty object
            % will be created
            
            if nargin == 1 && exist(filepath, 'file')
                metafile = load(filepath);
                obj.Tests = metafile.Tests;
            else
                obj.Tests = struct([]);
            end
        end
        
        function addRow(obj, test)
            %ADDROW Add a test            
            
            % Add missing fields to this object's Tests structure
            fields = setdiff(fieldnames(test), fieldnames(obj.Tests));
            emptyCell = cell(length(obj.Tests),1);
            for f = 1:length(fields)
                if any(strcmpi(fields{f}, fieldnames(obj.Tests)))
                    error(['Fieldnames must be unique and case-insensitive (%s)',...
                        fields{f}]);
                end
                [obj.Tests.(fields{f})] = emptyCell{:};
            end
            
            % Search for existing entries
            common = NeuroAnalysis.Base.getstructfields(test,...
                {'sourceformat', 'Subject_ID', 'RecordSession', ...
                'RecordSite', 'ID'});
            keys = fieldnames(common);
            values = struct2cell(common);
            roughMatch = obj.queryImpl(keys, values);
            if any(roughMatch)
                match = NeuroAnalysis.Base.EvalFun(...
                    ['NeuroAnalysis.',test.sourceformat,'.MatchMetadata'],...
                    {test, obj.Tests(roughMatch)});
            else
                match =  [];
            end
            
            % Merge if necessary
            if match
                % Use index of match to find index in roughMatch
                toMerge = find(roughMatch, match, 'first');
                toMerge = toMerge(end);
                
                % Merge the two tests
                newTest = NeuroAnalysis.Base.EvalFun(...
                    ['NeuroAnalysis.',test.sourceformat,'.MergeMetadata'],...
                    {obj.Tests(toMerge), test});
                
                % Remove the old row
                toKeep = ones(1,length(obj.Tests)); 
                toKeep(toMerge) = 0;
                obj.Tests = obj.Tests(logical(toKeep));
            else
                newTest = test;
            end
            
            % Add missing fields to the test
            fields = setdiff(fieldnames(obj.Tests), fieldnames(newTest));
            for f = 1:length(fields)
                newTest.(fields{f}) = {};
            end
            
            % Ready to merge
            % assert(isempty(setdiff(fieldnames(obj.Tests), fieldnames(newTest))));
            obj.Tests(end+1) = newTest;
        end
        
        function mt = query(obj, varargin)
            %QUERY Find matching tests
            nargs = length(varargin);
            if nargs == 0 || mod(nargs, 2) == 1
               error(['inputs should be name,value pairs, e.g. ',...
                   'query(MetaTable, ''Subject_ID'', ''R1701'')'])
            end
            keys = varargin(1:2:end);
            values = varargin(2:2:end);
            match = obj.queryImpl(keys, values);
            mt = NeuroAnalysis.IO.MetaTable();
            mt.Tests = obj.Tests(match);
        end
        
        function Tests = list(obj)
            %LIST list all the tests
            Tests = obj.Tests;
        end
        
        function export(obj, filepath)
            %EXPORT save to disk
            disp(['Saving metadata:    ',filepath,'    ...']);
            metafile.Tests = obj.Tests;
            fields = fieldnames(obj.Tests);
            for f=1:length(fields)
                metafile.TestsScalar.(fields{f}) = {obj.Tests.(fields{f})};
            end
            save(filepath, '-struct', 'metafile', '-v7.3');
            disp('Saving metadata:    Done.');
        end
        
        function [Invalid] = synchronize(obj, verbose)
            %SYNCHRONIZE remove rows with missing files
            % Returns a struct array of invalid tests
            disp('Synchronizing metadata:   ...');
            Invalid = struct([]);
            fields = [fieldnames(obj.Tests); {'missingFiles'}];
            emptyCell = cell(0,1);
            for f = 1:length(fields)
                [Invalid.(fields{f})] = emptyCell{:};
            end
            if isempty(obj.Tests)
                disp('Synchronizing metadata:   Done. (Empty metadata file)');
                return;
            end
            filesExist = true(length(obj.Tests),1); 
            for t = 1:length(obj.Tests)
                test = obj.Tests(t);
                test.missingFiles = {};
                % Check that some files are present
                if isempty(test.files)
                    filesExist(t) = false;
                end
                % Check that all files exist
                for f = 1:length(test.files)
                    if ~exist(test.files{f},'file')
                        filesExist(t) = false;
                        test.missingFiles = [test.missingFiles test.files(f)];
                        if verbose
                            disp(['Synchronizing metadata:    Missing file: ',...
                                test.files{f}]);
                        end     
                    end
                end
                Invalid(t) = test;
            end
            Invalid = Invalid(~filesExist);
            obj.Tests = obj.Tests(filesExist);
            disp('Synchronizing metadata:    Done.');
        end
    end
    
    methods (Access = private)
        
        function match = queryImpl(obj, keys, values)
            %QUERYIMPL Internal query for matching tests
            match = true(length(obj.Tests),1);
            for t = 1:length(obj.Tests)
                for kv = 1:length(keys)
                    match(t) = match(t) && ...
                        isequal(obj.Tests(t).(keys{kv}), values{kv});
                end
            end
        end
        
    end
    
end

