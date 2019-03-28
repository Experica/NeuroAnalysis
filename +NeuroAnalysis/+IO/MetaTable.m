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
            
            % Add UUID for test
            test.UUID = char(java.util.UUID.randomUUID());
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
            if isempty(searchtemplate)
                keys = {};
                values = {};
            else
                keys = fieldnames(searchtemplate);
                values = struct2cell(searchtemplate);
            end
            roughmatchindex = obj.iquery(keys, values);
            if ~isempty(roughmatchindex)
                matchindex = NeuroAnalysis.Base.EvalFun(...
                    ['NeuroAnalysis.',test.sourceformat,'.FindMetadata'],...
                    {obj,test,roughmatchindex});
            end
            
            if ~isempty(matchindex)
                % Keep old UUID
                for i = 1:length(matchindex)
                    uuid = obj.Tests(matchindex(i)).UUID;
                    if ~isempty(uuid)
                        test.UUID=uuid;
                        break;
                    end
                end
                % Merge new test and duplicate
                test= NeuroAnalysis.Base.EvalFun(...
                    ['NeuroAnalysis.',test.sourceformat,'.MergeMetadata'],...
                    {obj, test,matchindex(1)});
                % Delete duplicate rows
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
            for ii = 1:length(range)
                i=range(ii);
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
            if dosort && ~isempty(index)
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
            disp(['Saving metadata:    ',filepath,'    Done.']);
        end
        
        function [missingtest] = validate(obj, dataroot,verbose)
            %VALIDATE ensure metatable contains valid data files relative to DataRoot directory.
            %Returns a struct array of tests where invalid data files appeared
            
            if nargin ==2
                verbose=false;
            end
            deletemissingfile = false;
            forall = false;
            
            missingtest=[];
            disp('Validating metadata:    ...');
            if isempty(obj.Tests)
                disp('Validating metadata:   Empty.');
                return;
            end
            
            % Check for missing files
            missingindex=[];
            for t = 1:length(obj.Tests)
                test = obj.Tests(t);
                missingfiles = {};
                
                if ~isempty(test.files)
                    for i = length(test.files):-1:1
                        if startsWith(test.files{i}, '.')
                            p=fullfile(dataroot,test.files{i});
                        else
                            p=test.files{i};
                        end
                        if ~exist(p,'file')
                            if verbose
                                disp(['Validating metadata:    Missing File: ',p]);
                            end
                            missingfiles = [missingfiles p];
                            
                            if ~forall
                                choice = questdlg(['Remove Missing File: ', p,' ?'],'Missing File Action',...
                                    'No(All)','No','Yes','No(All)');
                                switch choice
                                    case 'Yes'
                                        deletemissingfile = true;
                                    case 'Yes(All)'
                                        deletemissingfile = true;
                                        forall = true;
                                    case 'No'
                                        deletemissingfile = false;
                                    case 'No(All)'
                                        deletemissingfile = false;
                                        forall = true;
                                end
                            end
                            if deletemissingfile
                                obj.Tests(t).files(i)=[];
                            end
                        else
                            obj.Tests(t).files{i}= fullfile('.',strrep(p,dataroot,''));
                        end
                    end
                end
                if isempty(obj.Tests(t).files)
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
                        disp(['Validating metadata:    Empty Column: ',...
                            fields{f},' will be removed.']);
                    end
                    emptyfields = [emptyfields fields(f)];
                end
            end
            obj.Tests = rmfield(obj.Tests, emptyfields);
            
            disp('Validating metadata:    Done.');
        end
        
    end
    
end