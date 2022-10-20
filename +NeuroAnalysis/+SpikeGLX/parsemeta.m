function [meta] = parsemeta(meta,d,probeinfo)
%PREPAREMETA Summary of this function goes here
%   Detailed explanation goes here

if nargin == 1
    d='ap';
    probeinfo=[];
end

if strcmp(meta.typeThis,'imec')
    meta.fs = meta.imSampRate;
    meta.snsApLfSy = int64(meta.snsApLfSy);
    meta.acqApLfSy = int64(meta.acqApLfSy);
    % imec probe version
    if isfield(meta, 'imDatPrb_type')
        meta.probeversion = meta.imDatPrb_type;
    else
        meta.probeversion = -1; % Phase3A
    end
    switch (meta.probeversion)
        case {21,24}
            % factor for converting 16-bit data to voltage
            fi2v = meta.imAiRangeMax / 8192;
        otherwise % Phase3A - 1.0
            % probe channel spacing[x,y,z] in um
            spacing = [32,20,0];
            % factor for converting 16-bit data to voltage
            fi2v = meta.imAiRangeMax / 512;
            % No. of channels multiplexed into one ADC
            nchmx = 12;
            % ADC1 {0,2,4,6,8,10,12,14,16,18,20,22}
            % ADC2 {1,3,5,7,9,11,13,15,17,19,21,23}
            % ADC3 {24,26,28,30,32,34,36,38,40,42,44,46}
            % ADC4 {25,27,29,31,33,35,37,39,41,43,45,47}
            % ...
            nch = meta.acqApLfSy(1);
            dmxgroup = cell(nchmx,1);
            hc = 1:2*nchmx:nch;
            dmxgroup{1} = [hc,hc+1];
            for g = 2:nchmx
                dmxgroup{g} = dmxgroup{g-1}+2;
            end
            % reference IDs
            if isfield(meta, 'imProbeOpt') % Phase3A
                switch  meta.imProbeOpt
                    case 4
                        refch = int64([36, 75, 112, 151, 188, 227, 264]+1);
                    otherwise
                        refch = int64([36, 75, 112, 151, 188, 227, 264, 303, 340, 379]+1);
                end
                probesn = meta.imProbeSN;
                rofmt = '(%d %d %d %d %d';
            else % Phase3B - 1.0
                refch = int64(191+1); % 192, 576, 960 for bank 0, 1, 2
                probesn = meta.imDatPrb_sn;
                rofmt = '(%d %d %d %d %d %d';
            end
            % imec readout table
            C = textscan(meta.imroTbl, rofmt, ...
                'EndOfLine', ')', 'HeaderLines', 1 );
            meta.roch = int64(cell2mat(C(1)));
            meta.robank = int64(cell2mat(C(2)));
            meta.rorefch = int64(cell2mat(C(3)));
            meta.roapgain = cell2mat(C(4));
            meta.rolfgain = cell2mat(C(5));
    end
    meta.nchmx = nchmx;
    meta.dmxgroup = dmxgroup;
    meta.probespacing = spacing;
    meta.refch = refch;
    meta.fi2v = fi2v;
    meta.probesn = probesn;
    % bad channels
    psn = ['SN_',num2str(meta.probesn)];
    if isstruct(probeinfo) && isfield(probeinfo,psn)
        badch = probeinfo.(psn).(d).badch;
        for i = 1:length(badch)
            if ischar(badch{i})
                badch{i} = str2num(badch{i});
            end
        end
        meta.badch = int64(unique(cell2mat(badch))+1);
    end
    % exclude channels, here only handle external referencing
    if all(meta.rorefch==0)
        excludechans = meta.refch;
    end
    if isfield(meta,'badch')
        excludechans = union(excludechans,meta.badch);
    end
    meta.excludechans = excludechans;
    % imec shank map for saved channels
    if isfield(meta,'snsShankMap')
        header = int64(str2num(regexp(meta.snsShankMap,'([0-9,]*)','match','once')));
        meta.nshank = header(1);
        meta.ncol= header(2);
        meta.nrow = header(3);
        C = textscan(meta.snsShankMap, '(%d:%d:%d:%*s', ...
            'EndOfLine', ')', 'HeaderLines', 1 );
        meta.savedshanks = int64(cell2mat(C(1))+1);
        meta.savedcols = int64(cell2mat(C(2))+1);
        meta.savedrows = int64(cell2mat(C(3))+1);
        meta.nshanksaved = int64(length(unique(meta.savedshanks)));
        meta.ncolsaved = int64(length(unique(meta.savedcols)));
        meta.nrowsaved = int64(length(unique(meta.savedrows)));
    end
    meta.syncch = 6+1; % fixed bit 6 for imec sync in PXI system
else % nidq
    meta.fs = meta.niSampRate;
    % factor for converting 16-bit data to voltage
    meta.fi2v = meta.niAiRangeMax / 32768;
    meta.snsMnMaXaDw = int64(meta.snsMnMaXaDw);
    meta.syncch = meta.syncNiChan + 1;
end
% Return original channel IDs, because the ith channel in the file isn't necessarily
% the ith acquired channel, so it could be used to index ith saved to original.
if ischar(meta.snsSaveChanSubset) && strcmp(meta.snsSaveChanSubset, 'all')
    meta.savedchans = int64(1:meta.nSavedChans);
else
    meta.savedchans = int64(meta.snsSaveChanSubset+1);
end

end

