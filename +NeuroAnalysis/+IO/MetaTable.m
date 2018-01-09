classdef MetaTable < handle
    %METATABLE Metadata table
    %   Detailed explanation goes here
    
    properties (Access = private)
        Tests       % Structure of metadata
        fileID
    end
    
    methods (Access = public)
        
        function obj = MetaTable(filepath)
            %METATABLE Create a MetaTable object from the given path
            % If the path does not point to a file, then an empty object
            % will be created
            if exist(filepath, 'file')
                metafile = load(filepath);
                obj.fileID = fopen(filepath); % to prevent access
                obj.Tests = metafile.Tests;
            else
                obj.Tests = struct([]);
                metafile = struct();
                save(filepath, '-struct', 'metafile', '-v7.3');
                obj.fileID = fopen(filepath);
            end
        end
        
        function addEntry(obj, test)
            %ADDENTRY Add a test            

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
            sf = false(length(obj.Tests),1);
            for t = 1:length(obj.Tests)
                sf(t) = isequal(obj.Tests(t).sourceformat, test.sourceformat);
            end
            match = NeuroAnalysis.Base.EvalFun(...
                ['NeuroAnalysis.',test.sourceformat,'.MatchMetadata'],...
                {test, obj.Tests(sf)});
            
            % Merge if necessary
            if match
                toMerge = find(sf, match, 'first');
                toMerge = toMerge(end);
                newTest = NeuroAnalysis.Base.EvalFun(...
                    ['NeuroAnalysis.',test.sourceformat,'.MergeMetadata'],...
                    {test, obj.Tests(toMerge)});
                
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
            
            obj.Tests(end+1) = newTest;
        end
        
        function Tests = query(obj, varargin)
            %QUERY Find matching tests
            nargs = length(varargin);
            if nargs == 0 || mod(nargs, 2) == 1
               error('inputs should be name,value pairs, e.g. query(MetaTable, ''Subject_ID'', ''R1701'')')
            end
            keys = varargin(1:2:end);
            values = varargin(2:2:end);
            isMatch = true(length(obj.Tests),1);
            for t = 1:length(obj.Tests)
                for kv = 1:length(keys)
                    isMatch(t) = isMatch(t) && isequal(obj.Tests(t).(keys{kv}), values{kv});
                end
            end
            Tests = obj.Tests(isMatch);
        end
        
        function Tests = list(obj)
            %LIST list all the tests
            Tests = obj.Tests;
        end
        
        function export(obj, filepath)
            %EXPORT save to disk
            disp(['Saving metadata:    ',filepath,'    ...']);
            fclose(obj.fileID);
            metafile.Tests = obj.Tests;
            fields = fieldnames(obj.Tests);
            for f=1:length(fields)
                metafile.TestsScalar.(fields{f}) = {obj.Tests.(fields{f})};
            end
            save(filepath, '-struct', 'metafile', '-v7.3');
            disp('Saving metadata:    Done.');
        end
        
        function [Missing] = synchronize(obj, verbose)
            %SYNCHRONIZE remove entries without existing files
            disp('Synchronizing metadata:   ...');
            Missing = struct([]);
            fields = [fieldnames(obj.Tests); {'missingFiles'}];
            emptyCell = cell(0,1);
            for f = 1:length(fields)
                [Missing.(fields{f})] = emptyCell{:};
            end
            if isempty(obj.Tests)
                return;
            end
            filesExist = true(length(obj.Tests),1); 
            for t = 1:length(obj.Tests)
                test = obj.Tests(t);
                test.missingFiles = {};
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
                Missing(t) = test;
            end
            Missing = Missing(~filesExist);
            obj.Tests = obj.Tests(filesExist);
            disp('Synchronizing metadata:    Done.');
        end
    end
     
    methods (Static, Access = private)

        function value = getfieldi(S,field)
            %GETFIELDI case-insensitive dynamic structure field access
            names = fieldnames(S);
            isField = strcmpi(field,names);  

            if any(isField)
                assert(sum(isField) == 1)
                value = {S.(names{isField})};
            else
                value = {};
            end
        end
        
    end
    
end

