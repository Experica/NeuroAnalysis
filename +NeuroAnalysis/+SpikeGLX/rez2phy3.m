function rez2phy3(rez, savePath,dataset,d)
% pull out results from kilosort's rez to either return to workspace or to
% save in the appropriate format for the phy GUI to run on. If you provide
% a savePath it should be a folder, and you will need to have npy-matlab
% available (https://github.com/kwikteam/npy-matlab)
%
% This is a modified version of the `rezToPhy2` function in Kilosort3

[~, Nfilt, Nrank] = size(rez.W);
% for Phy, we need to pad the spikes with zeros so the spikes are aligned to the center of the window
rez.Wphy = cat(1, zeros(1+rez.ops.nt0min, Nfilt, Nrank), rez.W);

rez.W = gather(single(rez.Wphy));
rez.U = gather(single(rez.U));
rez.mu = gather(single(rez.mu));

if size(rez.st3,2)>4
    rez.st3 = rez.st3(:,1:4);
end

[~, isort]   = sort(rez.st3(:,1), 'ascend');
rez.st3      = rez.st3(isort, :);
if ~isempty(rez.cProj)
    rez.cProj    = rez.cProj(isort, :);
    rez.cProjPC  = rez.cProjPC(isort, :, :);
end

fs = dir(fullfile(savePath, '*.npy'));
for i = 1:length(fs)
    delete(fullfile(savePath, fs(i).name));
end
if exist(fullfile(savePath, '.phy'), 'dir')
    rmdir(fullfile(savePath, '.phy'), 's');
end

% spikeTimes will be in samples, not seconds
spikeTimes = uint64(rez.st3(:,1));
% account for ops.trange(1) to accomodate real time
spikeTimes = spikeTimes - rez.ops.trange(1)*rez.ops.fs;
spikeTemplates = uint32(rez.st3(:,2));
if size(rez.st3,2)>4
    spikeClusters = uint32(1+rez.st3(:,5));
end
amplitudes = rez.st3(:,3);

Nchan = rez.ops.Nchan;
xcoords     = rez.xcoords(:);
ycoords     = rez.ycoords(:);
chanMap     = rez.ops.chanMap(:);
chanMap0ind = chanMap - 1;
U = rez.U;
W = rez.W;
nt0 = size(W,1);
Nfilt = size(W,2);

templates = zeros(Nchan, nt0, Nfilt, 'single');
for iNN = 1:Nfilt
    templates(:,:,iNN) = squeeze(U(:,iNN,:)) * squeeze(W(:,iNN,:))';
end
templates = permute(templates, [3 2 1]); % now it's nTemplates x nSamples x nChannels
templatesInds = repmat([0:size(templates,3)-1], size(templates,1), 1); % we include all channels so this is trivial

templateFeatures = rez.cProj;
templateFeatureInds = uint32(rez.iNeigh);
pcFeatures = rez.cProjPC;
pcFeatureInds = uint32(rez.iNeighPC);

whiteningMatrix_raw = rez.Wrot/rez.ops.scaleproc;
whiteningMatrixInv_raw = whiteningMatrix_raw^-1;
whiteningMatrix = eye(size(rez.Wrot)) / rez.ops.scaleproc;
whiteningMatrixInv = whiteningMatrix^-1;

% unwhiten all the templates
tempsUnW = zeros(size(templates));
for t = 1:size(templates,1)
    tempsUnW(t,:,:) = squeeze(templates(t,:,:))*whiteningMatrixInv;
end

% The amplitude on each channel is the positive peak minus the negative
tempChanAmps = squeeze(max(tempsUnW,[],2))-squeeze(min(tempsUnW,[],2));

% The template amplitude is the amplitude of its largest channel
tempAmpsUnscaled = max(tempChanAmps,[],2);

% assign all spikes the amplitude of their template multiplied by their
% scaling amplitudes
spikeAmps = tempAmpsUnscaled(spikeTemplates).*amplitudes;

% take the average of all spike amps to get actual template amps (since
% tempScalingAmps are equal mean for all templates)
ta = clusterAverage(spikeTemplates, spikeAmps);
tids = unique(spikeTemplates);
tempAmps = zeros(numel(rez.mu),1);
tempAmps(tids) = ta; % because ta only has entries for templates that had at least one spike
gain = getOr(rez.ops, 'gain', 1);
tempAmps = gain*tempAmps'; % for consistency, make first dimension template number

if ~isempty(savePath)
    writeNPY(spikeTimes, fullfile(savePath, 'spike_times.npy'));
    writeNPY(uint32(spikeTemplates-1), fullfile(savePath, 'spike_templates.npy')); % -1 for zero indexing
    if size(rez.st3,2)>4
        writeNPY(uint32(spikeClusters-1), fullfile(savePath, 'spike_clusters.npy')); % -1 for zero indexing
    else
        writeNPY(uint32(spikeTemplates-1), fullfile(savePath, 'spike_clusters.npy')); % -1 for zero indexing
    end
    writeNPY(amplitudes, fullfile(savePath, 'amplitudes.npy'));
    writeNPY(templates, fullfile(savePath, 'templates.npy'));
    writeNPY(templatesInds, fullfile(savePath, 'templates_ind.npy'));
    
    writeNPY(int32(chanMap0ind), fullfile(savePath, 'channel_map_raw.npy'));
    writeNPY(int32(0:rez.ops.Nchan-1), fullfile(savePath, 'channel_map.npy'));
    writeNPY([xcoords ycoords], fullfile(savePath, 'channel_positions.npy'));
    
    if ~isempty(templateFeatures)
        writeNPY(templateFeatures, fullfile(savePath, 'template_features.npy'));
        writeNPY(templateFeatureInds'-1, fullfile(savePath, 'template_feature_ind.npy'));% -1 for zero indexing
        writeNPY(pcFeatures, fullfile(savePath, 'pc_features.npy'));
        writeNPY(pcFeatureInds'-1, fullfile(savePath, 'pc_feature_ind.npy'));% -1 for zero indexing
    end
    
    writeNPY(whiteningMatrix, fullfile(savePath, 'whitening_mat.npy'));
    writeNPY(whiteningMatrixInv, fullfile(savePath, 'whitening_mat_inv.npy'));
    writeNPY(whiteningMatrix_raw, fullfile(savePath, 'whitening_mat_raw.npy'));
    writeNPY(whiteningMatrixInv_raw, fullfile(savePath, 'whitening_mat_inv_raw.npy'));
    
    if isfield(rez, 'simScore')
        similarTemplates = rez.simScore;
        writeNPY(similarTemplates, fullfile(savePath, 'similar_templates.npy'));
    end
    
    % save a list of "good" clusters for Phy
    fileID = fopen(fullfile(savePath, 'cluster_KSLabel.tsv'),'w');
    fprintf(fileID, 'cluster_id%sKSLabel', char(9));
    fprintf(fileID, char([13 10]));
    
    fileIDCP = fopen(fullfile(savePath, 'cluster_ContamPct.tsv'),'w');
    fprintf(fileIDCP, 'cluster_id%sContamPct', char(9));
    fprintf(fileIDCP, char([13 10]));
    
    fileIDA = fopen(fullfile(savePath, 'cluster_Amplitude.tsv'),'w');
    fprintf(fileIDA, 'cluster_id%sAmplitude', char(9));
    fprintf(fileIDA, char([13 10]));
    
    rez.est_contam_rate(isnan(rez.est_contam_rate)) = 1;
    for j = 1:length(rez.good)
        if rez.good(j)
            fprintf(fileID, '%d%sgood', j-1, char(9));
        else
            fprintf(fileID, '%d%smua', j-1, char(9));
        end
        fprintf(fileID, char([13 10]));
        
        fprintf(fileIDCP, '%d%s%.1f', j-1, char(9), rez.est_contam_rate(j)*100);
        fprintf(fileIDCP, char([13 10]));
        
        fprintf(fileIDA, '%d%s%.1f', j-1, char(9), tempAmps(j));
        fprintf(fileIDA, char([13 10]));
    end
    fclose(fileID);
    fclose(fileIDCP);
    fclose(fileIDA);
    
    % Duplicate "KSLabel" as "group", a special metadata ID for Phy, so that
    % filtering works as expected in the cluster view
    KSLabelFilename = fullfile(savePath, 'cluster_KSLabel.tsv');
    copyfile(KSLabelFilename, fullfile(savePath, 'cluster_group.tsv'));
    
    % make params file
    if ~exist(fullfile(savePath,'params.py'),'file')
        fid = fopen(fullfile(savePath,'params.py'), 'w');
        % save dataset path so that phy sorting result could be loaded back into dataset
        if iscell(dataset.filepath)
            datasetpath = strjoin(dataset.filepath,', ');
        else
            datasetpath = dataset.filepath;
        end
        fprintf(fid,['dataset_path = ''',replace(datasetpath,'\','\\'),'''\n']);
        fprintf(fid,'secondperunit = %.32f\n',dataset.secondperunit);
        if isfield(dataset,'binfilerange')
            writeNPY(dataset.binfilerange, fullfile(savePath, 'binfilerange.npy'));
        end
        if isfield(rez,'dshift')
            writeNPY(rez.dshift, fullfile(savePath, 'dshift.npy'));
            writeNPY(rez.yblk, fullfile(savePath, 'yblk.npy'));
        end
        fprintf(fid,['rawdat_path = ''', strrep(rez.ops.fbinary, '\', '/') '''\n']);
        fprintf(fid,'n_channels_rawdat = %i\n',rez.ops.NchanTOT);
        fprintf(fid,['dat_path = ''', strrep(rez.ops.fproc, '\', '/') '''\n']);
        fprintf(fid,'n_channels_dat = %i\n',rez.ops.Nchan);
        fprintf(fid,'nsample = %i\n',dataset.(d).meta.nFileSamp);
        fprintf(fid,['imecindex = ''',d(3:end),'''\n']);
        fprintf(fid,'dtype = ''int16''\n');
        fprintf(fid,'offset = 0\n');
        fprintf(fid,'sample_rate = %.32f\n',rez.ops.fs);
        fprintf(fid,'sort_from = ''kilosort3''\n');
        if isfield(rez.ops,'fshigh')
            hp='True';
        else
            hp='False';
        end
        fprintf(fid,['hp_filtered = ', hp]);
        fclose(fid);
    end
    % save figures
    for n=1:length(findobj('type','figure'))
        f = figure(n);
        exportgraphics(f,fullfile(savePath,[num2str(n),'.png']));
        close(f);
    end
end

end

