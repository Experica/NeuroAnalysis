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
                {'Subject_ID', 'RecordSession'});
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
            matchindex = obj.iquery(keys, values, [], true);
            if isempty(matchindex)
                warning('No Matching Metadata Found. ');
                return;
            end
            mt = NeuroAnalysis.IO.MetaTable();
            mt.Tests = obj.Tests(matchindex);
        end
        
        function index = iquery(obj, keys, values, range, dosort)
            %IQUERY query for matching tests, return indices
            
           
            if nargin < 4 || isempty(range)
                range=1:length(obj.Tests);
            else
                range(range<1 & range>length(obj.Tests))=[];
            end
            if nargin < 5
                dosort=false;
            end
            
            index = [];
            dates = [];
            for i = range
                ismatch = true;
                for j = 1:length(keys)
                    ismatch = ismatch & isequal(obj.Tests(i).(keys{j}), values{j});
                end
                if ismatch
                    index=[index;i];
                    if dosort && isfield(obj.Tests, 'date') && ...
                            ~isempty(obj.Tests(i).date)
                        dates=[dates;obj.Tests(i).date];
                    else
                        dates=[dates;NaN];
                    end
                end
            end
            
            % Sort by date
            if dosort
                sorted = sortrows([dates, index],'descend',...
                    'MissingPlacement', 'last');
                index = sorted(:,2);
            end
            
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
            
            disp('Validating metadata:   ...');
            if isempty(obj.Tests)
                disp('Validating metadata:   Empty.');
                return;
            end
            
            % Check for missing files
            missingtest=[];
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
            
            % Check for empty fields
            fields = fieldnames(obj.Tests);
            emptyfields = {};
            for f = 1:length(fields)
                if all(cellfun(@isempty, {obj.Tests.(fields{f})}))
                    if verbose
                        disp(['Validating metadata:    Empty column: ',...
                            fields{f},' is removed.']);
                    end
                    emptyfields = [emptyfields fields(f)];
                end
            end
            obj.Tests = rmfield(obj.Tests, emptyfields);
            
            disp('Validating metadata:    Done.');
        end
        
    end
    
    methods (Access = private)
        
    end
    
end
