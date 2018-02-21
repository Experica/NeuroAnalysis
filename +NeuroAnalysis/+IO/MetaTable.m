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
                obj.Tests = [];
            end
        end
        
        function addRow(obj, test)
            %ADDROW Add a test
            
            % Add new fields to this object's Tests struct array
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
                obj.Test(matchindex)=[];
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
            %             fields = fieldnames(obj.Tests);
            %             for f=1:length(fields)
            %                 metafile.TestsScalar.(fields{f}) = {obj.Tests.(fields{f})};
            %             end
            save(filepath, '-struct', 'metafile', '-v7.3');
            disp('Saving metadata:    Done.');
        end
        
        function [missingtest] = synchronize(obj, verbose)
            %SYNCHRONIZE remove rows containing missing data files
            % Returns a struct array of invalid tests
            
            missingtest=[];
            disp('Synchronizing metadata:   ...');
            if isempty(obj.Tests)
                disp('Synchronizing metadata:   Done. (Empty Metadata)');
                return;
            end
            
            missingindex=[];
            for t = 1:length(obj.Tests)
                test = obj.Tests(t);
                missingFiles = {};
                
                if ~isempty(test.files)
                    for i = 1:length(test.files)
                        if ~exist(test.files{i},'file')
                            if verbose
                                disp(['Synchronizing metadata:    Missing file: ',...
                                    test.files{i},' is removed.']);
                            end
                            missingFiles = [missingFiles test.files(i)];
                            test.files(i)=[];
                            obj.Tests(t).files(i)=[];
                        end
                    end
                end
                if isempty(test.files)
                    missingindex=[missingindex;t];
                    test.missingFiles= missingFiles;
                    missingtest=[missingtest;test];
                elseif ~isempty(missingFiles)
                    test.missingFiles= missingFiles;
                    missingtest=[missingtest;test];
                end
            end
            obj.Tests(missingindex)=[];
            disp('Synchronizing metadata:    Done.');
        end
    end
    
    methods (Access = private)
        
        function index = iquery(obj, keys, values)
            %IQUERY internal query for matching tests, return indices
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

