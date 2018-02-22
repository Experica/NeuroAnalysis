classdef MetaTable < handle
    %METATABLE Metadata table
    %   Detailed explanation goes here
    
    properties (GetAccess = public, SetAccess = private)
        Tests       % metadata table in struct array
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
            
            % Add new fields to Tests struct array
            newfields = setdiff(fieldnames(test), fieldnames(obj.Tests));
            emptyColumn = cell(length(obj.Tests),1);
            for i = 1:length(newfields)
                [obj.Tests.(newfields{i})] = emptyColumn{:};
            end
            
            % Search for duplicate old rows
            matchindex =  [];
            searchtemplate = NeuroAnalysis.Base.getstructfields(test,...
                {'sourceformat', 'Subject_ID', 'RecordSession', ...
                'RecordSite', 'ID'});
            keys = fieldnames(searchtemplate);
            values = struct2cell(searchtemplate);
            roughmatchindex = obj.iquery(keys, values);
            if ~isempty(roughmatchindex)
                matchindex = NeuroAnalysis.Base.EvalFun(...
                    ['NeuroAnalysis.',test.sourceformat,'.FindMetadata'],...
                    {obj,test,roughmatchindex});
            end
            
            if ~isempty(matchindex)
                % Merge old row and test
                test= NeuroAnalysis.Base.EvalFun(...
                    ['NeuroAnalysis.',test.sourceformat,'.MergeMetadata'],...
                    {obj, test,matchindex});
                % Delete old row
                obj.Tests(matchindex)=[];
            else
                % Add missing fields to the test
                missingfields = setdiff(fieldnames(obj.Tests), fieldnames(test));
                for i = 1:length(missingfields)
                    test.(missingfields{i}) = {};
                end
            end
            % Add new row
            obj.Tests(end+1) = test;
        end
        
        function mt = query(obj, varargin)
            %QUERY Find matching metatable
            
            mt=[];
            nargs = length(varargin);
            if nargs == 0 || mod(nargs, 2) == 1
                warning(['inputs should be name,value pairs, e.g. ',...
                    'query(MetaTable, ''Subject_ID'', ''R1701'')'])
                return;
            end
            keys = varargin(1:2:end);
            values = varargin(2:2:end);
            matchindex = obj.iquery(keys, values);
            if isempty(matchindex)
                warning('No Matching Metadata Found. ');
                return;
            end
            mt = NeuroAnalysis.IO.MetaTable();
            mt.Tests = obj.Tests(matchindex);
        end
        
        function export(obj, filepath)
            %EXPORT save metadata
            
            disp(['Saving metadata:    ',filepath,'    ...']);
            metafile.Tests = obj.Tests;
            save(filepath, '-struct', 'metafile', '-v7.3');
            disp('Saving metadata:    Done.');
        end
        
        function [missingtest] = validate(obj, verbose)
            %VALIDATE ensure metatable contains valid data files
            % Returns a struct array of tests where invalid data files appeared
            
            if nargin ==1
                verbose=false;
            end
            
            missingtest=[];
            disp('Validating metadata:   ...');
            if isempty(obj.Tests)
                disp('Validating metadata:   Empty.');
                return;
            end
            
            missingindex=[];
            for t = 1:length(obj.Tests)
                test = obj.Tests(t);
                missingfiles = {};
                
                if ~isempty(test.files)
                    for i = 1:length(test.files)
                        if ~exist(test.files{i},'file')
                            if verbose
                                disp(['Validating metadata:    Missing file: ',...
                                    test.files{i},' is removed.']);
                            end
                            missingfiles = [missingfiles test.files(i)];
                            test.files(i)=[];
                            obj.Tests(t).files(i)=[];
                        end
                    end
                end
                if isempty(test.files)
                    missingindex=[missingindex;t];
                    test.missingfiles= missingfiles;
                    missingtest=[missingtest;test];
                elseif ~isempty(missingfiles)
                    test.missingfiles= missingfiles;
                    missingtest=[missingtest;test];
                end
            end
            obj.Tests(missingindex)=[];
            disp('Validating metadata:    Done.');
        end
    end
    
    methods (Access = private)
        
        function index = iquery(obj, keys, values)
            %IQUERY Internal query for matching tests, return indices
            
            index = [];
            for i = 1:length(obj.Tests)
                ismatch = true;
                for j = 1:length(keys)
                    ismatch = ismatch & isequal(obj.Tests(i).(keys{j}), values{j});
                end
                if ismatch
                    index=[index,i];
                end
            end
        end
        
    end
    
end

