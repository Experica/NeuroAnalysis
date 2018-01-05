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
            if exist(filepath, 'file')
                metafile = load(filepath);
                obj.Tests = metafile.Tests;
            else
                obj.Tests = struct([]);
            end
        end
        
        function addEntries(obj, tests, deleteMissing)
            %ADDENTRIES Add multiple tests
            % if deleteMissing is true, then any tests not supplied in the
            % input will be removed
            oldIdx = length(obj.Tests);
            for t = 1:length(tests)
                [obj.Tests, merged] = obj.merge(obj.Tests, tests{t});
                oldIdx = oldIdx - merged;
            end
            if deleteMissing
                obj.Tests = obj.Tests(oldIdx+1:end);
            end
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
            metafile.Tests = obj.Tests;
            fields = fieldnames(obj.Tests);
            for f=1:length(fields)
                metafile.TestsScalar.(fields{f}) = {obj.Tests.(fields{f})};
            end
            save(filepath, '-struct', 'metafile', '-v7.3');
            disp('Saving metadata:    Done.');
        end
    end
     
    methods (Static, Access = private)

         function [Tests, merged] = merge(Tests, test)
            %MERGE Merge a test to a structure of tests

            % Add missing fields to this object's Tests structure
            fields = setdiff(fieldnames(test), fieldnames(Tests));
            emptyCell = cell(length(Tests),1);
            for f = 1:length(fields)
                if any(strcmpi(fields{f}, fieldnames(Tests)))
                    error(['Fieldnames must be unique and case-insensitive (%s)',...
                        fields{f}]);
                end
                [Tests.(fields{f})] = emptyCell{:};
            end
            
            % Merge according to sourceformat if necessary
            toMerge = false(length(Tests),1);
            for t = 1:length(Tests)
                toMerge(t) = isequal(Tests(t).(test.key), test.(test.key));
            end
            merged = sum(toMerge);
            if merged
                assert(merged == 1)
                newTest = NeuroAnalysis.Base.EvalFun(...
                    ['NeuroAnalysis.',test.sourceformat,'.MetadataMerge'],...
                    {test, Tests(toMerge)});
            else
                newTest = test;
            end
            
            % Add missing fields to the test
            Tests = Tests(~toMerge);
            fields = setdiff(fieldnames(Tests), fieldnames(newTest));
            for f = 1:length(fields)
                newTest.(fields{f}) = {};
            end
            Tests(end+1) = newTest;
        end
        
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

