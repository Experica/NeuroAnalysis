function [ dataset ] = Prepare( filepath,varargin )
%PREPARE Read Ripple data by Ripple neuroshare API, VLab data and prepare dataset
%   Detailed explanation goes here

p = inputParser;
addRequired(p,'filepath');
addOptional(p,'datatype',{'Spike','LFP','Hi-Res','Stim','Analog30k','Analog1k','Digital'}); %'Raw'
addOptional(p,'electroderange',1:5120);
addOptional(p,'analogrange',10241:10270);
parse(p,filepath,varargin{:});
filepath = p.Results.filepath;
datatype = p.Results.datatype;
electroderange = p.Results.electroderange;
analogrange = p.Results.analogrange;

import NeuroAnalysis.Ripple.*
%% Prepare all data files
dataset = struct([]);
[hfile] = fopen(filepath,'r');
if hfile == -1
    warn(['Can not open file: ',filepath]);
    return;
end
[basepath,filename,ext] = fileparts(filepath);

isrippledata= false;
datafilepath = fullfile(basepath,filename);
[ns_RESULT, hFile] = ns_OpenFile(datafilepath);
if(strcmp(ns_RESULT,'ns_OK'))
    isrippledata = true;
end

isvlabdata = false;
vlabfilepath =fullfile(basepath,[filename '.yaml']);
if(exist(vlabfilepath,'file')==2)
    isvlabdata = true;
end

isvisstimdata = false;
visstimfilepath =fullfile(basepath,[filename '.mat']);
if(exist(visstimfilepath,'file')==2)
    isvisstimdata = true;
end
%% Read Ripple data
disp(['Reading Ripple Files:    ',datafilepath,'.*    ...']);
if(isrippledata)
    [ns_RESULT, nsFileInfo] = ns_GetFileInfo(hFile);
    if(~strcmp(ns_RESULT,'ns_OK'))
        ns_RESULT = ns_CloseFile(hFile);
        return;
    end
    
    EntityID = struct;
    EntityFileType = arrayfun(@(e)e.FileType,hFile.Entity);
    EntityType = arrayfun(@(e)e.EntityType,hFile.Entity,'Uniformoutput',false);
    EntityReason = arrayfun(@(e)e.Reason,hFile.Entity,'Uniformoutput',false);
    EntityElectrodeID = arrayfun(@(e)e.ElectrodeID,hFile.Entity);
    for i = 1:length(hFile.FileInfo)
        switch hFile.FileInfo(i).Type
            case 'nev'
                if ismember('Spike',datatype)
                    entityid = find(EntityFileType==i);
                    electrodeid = EntityElectrodeID(entityid);
                    vch = ismember(electrodeid,electroderange);
                    EntityID.Spike = entityid(vch);
                    ElectrodeID.Spike = electrodeid(vch);
                end
                if ismember('Digital',datatype)
                    EntityID.Digital = find((EntityFileType==i)&(cellfun( @(x)strcmp(x,'Event'),EntityType)));
                    Reason = EntityReason(EntityID.Digital);
                end
            case 'ns2'
                if ismember('LFP',datatype)
                    entityid = find(EntityFileType==i);
                    electrodeid = EntityElectrodeID(entityid);
                    vch = ismember(electrodeid,electroderange);
                    if any(vch)
                        EntityID.LFP = entityid(vch);
                        ElectrodeID.LFP = electrodeid(vch);
                    end
                    ns2TimeStamps = hFile.FileInfo(i).TimeStamps;
                end
                if ismember('Analog1k',datatype)
                    entityid = find(EntityFileType==i);
                    electrodeid = EntityElectrodeID(entityid);
                    vch = ismember(electrodeid,analogrange);
                    if any(vch)
                        EntityID.Analog1k = entityid(vch);
                        ElectrodeID.Analog1k = electrodeid(vch);
                    end
                    ns2TimeStamps = hFile.FileInfo(i).TimeStamps;
                end
            case 'ns5'
                if ismember('Raw',datatype)
                    entityid = find(EntityFileType==i);
                    electrodeid = EntityElectrodeID(entityid);
                    vch = ismember(electrodeid,electroderange);
                    if any(vch)
                        EntityID.Raw = entityid(vch);
                        ElectrodeID.Raw = electrodeid(vch);
                    end
                    ns5TimeStamps = hFile.FileInfo(i).TimeStamps;
                end
                if ismember('Analog30k',datatype)
                    entityid = find(EntityFileType==i);
                    electrodeid = EntityElectrodeID(entityid);
                    vch = ismember(electrodeid,analogrange);
                    if any(vch)
                        EntityID.Analog30k = entityid(vch);
                        ElectrodeID.Analog30k = electrodeid(vch);
                    end
                    ns5TimeStamps = hFile.FileInfo(i).TimeStamps;
                end
        end
    end
    
    dataset=struct;
    fdatatype = fieldnames(EntityID);
    for f=1:length(fdatatype)
        switch fdatatype{f}
            case 'Spike'
                for e=1:length(ElectrodeID.Spike)
                    [ns_RESULT, nsEntityInfo] = ns_GetEntityInfo(hFile, EntityID.Spike(e));
                    
                    spike = struct;
                    for i = 1:nsEntityInfo.ItemCount
                        [ns_RESULT, spike.time(i), spike.data(:,i), ~, spike.unitid(i)] = ns_GetSegmentData(hFile, EntityID.Spike(e), i);
                    end
                    spike.electrodeid = ElectrodeID.Spike(e);
                    dataset.spike(e) = spike;
                end
            case 'Digital'
                for e=1:length(Reason)
                    [ns_RESULT, nsEntityInfo] = ns_GetEntityInfo(hFile, EntityID.Digital(e));
                    
                    digital = struct;
                    for i = 1:nsEntityInfo.ItemCount
                        [ns_RESULT, digital.time(i), digital.data(i), ~] = ns_GetEventData(hFile, EntityID.Digital(e), i);
                    end
                    digital.channel = Reason(e);
                    dataset.digital(e) = digital;
                end
            case 'LFP'
                [ns_RESULT, Data] = ns_GetAnalogDataBlock(hFile, EntityID.LFP, 1, ns2TimeStamps(end)-ns2TimeStamps(1));
                [ns_RESULT, nsAnalogInfo] = ns_GetAnalogInfo(hFile, EntityID.LFP(1));
                
                dataset.lfp.data = Data;
                dataset.lfp.fs = nsAnalogInfo.SampleRate;
                dataset.lfp.electrodeid = ElectrodeID.LFP;
                dataset.lfp.time = (ns2TimeStamps/nsAnalogInfo.SampleRate);
            case 'Analog1k'
                [ns_RESULT, Data] = ns_GetAnalogDataBlock(hFile, EntityID.Analog1k, 1, ns2TimeStamps(end)-ns2TimeStamps(1));
                [ns_RESULT, nsAnalogInfo] = ns_GetAnalogInfo(hFile, EntityID.Analog1k(1));
                
                dataset.analog1k.data = Data;
                dataset.analog1k.fs = nsAnalogInfo.SampleRate;
                dataset.analog1k.electrodeid = ElectrodeID.Analog1k;
                dataset.analog1k.time = (ns2TimeStamps/nsAnalogInfo.SampleRate);
            case 'Raw'
                [ns_RESULT, Data] = ns_GetAnalogDataBlock(hFile, EntityID.Raw, 1, ns5TimeStamps(end)-ns5TimeStamps(1));
                [ns_RESULT, nsAnalogInfo] = ns_GetAnalogInfo(hFile, EntityID.Raw(1));
                
                dataset.raw.data = Data;
                dataset.raw.fs = nsAnalogInfo.SampleRate;
                dataset.raw.electrodeid = ElectrodeID.Raw;
                dataset.raw.time = (ns5TimeStamps/nsAnalogInfo.SampleRate);
            case 'Analog30k'
                [ns_RESULT, Data] = ns_GetAnalogDataBlock(hFile, EntityID.Analog30k, 1, ns5TimeStamps(end)-ns5TimeStamps(1));
                [ns_RESULT, nsAnalogInfo] = ns_GetAnalogInfo(hFile, EntityID.Analog30k(1));
                
                dataset.analog30k.data = Data;
                dataset.analog30k.fs = nsAnalogInfo.SampleRate;
                dataset.analog30k.electrodeid = ElectrodeID.Analog30k;
                dataset.analog30k.time = (ns5TimeStamps/nsAnalogInfo.SampleRate);
        end
    end
end

if ~isempty(dataset)
    dataset.source = datafilepath;
    dataset.sourceformat = 'Ripple';
end
ns_RESULT = ns_CloseFile(hFile);
disp('Reading Ripple Files:    Done.');
%% Prepare Ripple data

%% Prepare experimental data
if(isvlabdata && ~isempty(dataset))
    vlabdataset = NeuroAnalysis.VLab.Prepare(vlabfilepath);
    if ~isempty(vlabdataset)
        dataset.ex = vlabdataset.ex;
    end
elseif(isvisstimdata && ~isempty(dataset))
    visstimdataset = NeuroAnalysis.VisStim.Prepare(visstimfilepath);
    if ~isempty(visstimdataset)
        dataset.ex = visstimdataset.ex;
    end
end

end