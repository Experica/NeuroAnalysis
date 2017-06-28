function [ result ] = PlotData( exportpath, filename, sourceformat, dataset )
%PlotData Forward the dataset to tuning curve plugin
%   Detailed explanation goes here

switch sourceformat
    case 'Ripple'
        
        % Parameters
        whichElectrodes = [];
        dataPath = fileparts(exportpath);
        figuresPath = dataPath;
        fileName = filename;
        summaryFig = 0;
        plotLFP = 0;
        showFigures = 0;
        
        % Load what's needed
        fprintf('Loading data for plotting:\t')
        Params = dataset.ex.Params;
        [ Electrodes, LFP, AnalogIn, Events ] = convertDataset( dataset );
        fprintf('Done\n')
        
        % Adjust StimTimes
        pathparts = split(dataPath, filesep);
        Params.unitNo = sscanf(char(pathparts(end)), 'Unit%d');
        [stimOnTimes, stimOffTimes, source, latency, variation, hasError, errorMsg] = ...
            adjustStimTimes(Params, Events);
        StimTimes = [];
        StimTimes.on = stimOnTimes;
        StimTimes.off = stimOffTimes;
        StimTimes.latency = latency;
        StimTimes.variation = variation;
        StimTimes.source = source;
        
        if hasError
            plotStimTimes(StimTimes, AnalogIn, figuresPath, fileName, errorMsg);
        end
        
        if hasError > 1
            warning('Loading experiment failed. Skipping.');
            return;
        end
        
        % Number of bins or size of bins?
        switch Params.stimType
            case {'VelocityConstantCycles', 'VelocityConstantCyclesBar',...
                    'Looming', 'ConstantCycles'}
                Params.binType = 'number'; % constant number of bins
            otherwise
                Params.binType = 'size'; % constant size of bins (for F1)
        end
        
        % What to plot?
        switch Params.stimType
            case {'LatencyTest', 'LaserON', 'LaserGratings', ...
                    'NaturalImages', 'NaturalVideos'}
                plotMaps = 0;
                plotTCs = 0;
                plotBars = 1;
                plotRasters = 1;
                plotWFs = 1;
                plotISI = 1;
            case {'RFmap', 'CatRFdetailed', 'CatRFfast', 'CatRFfast10x10'}
                plotMaps = 1;
                plotTCs = 0;
                plotBars = 0;
                plotRasters = 1;
                plotWFs = 0;
                plotISI = 0;
            otherwise
                plotMaps = 0;
                plotTCs = 1;
                plotBars = 0;
                plotRasters = 1;
                plotWFs = 0;
                plotISI = 0;
        end
        
        
        % Do the analysis / plotting
        switch Params.stimType
            case {'OriLowHighTwoApertures', 'CenterNearSurround'}
                fprintf(2, 'Center surround not yet implemented.\n');
            otherwise
                
                % Binning and making rastergrams
                disp('analyzing...');
                [Params, Results] = analyze(dataPath, fileName, ...
                    Params, StimTimes, Electrodes, whichElectrodes);
                
                % Plotting tuning curves and maps
                disp('plotting...');
                plotAllResults(figuresPath, fileName, Params, Results, ...
                    whichElectrodes, plotTCs, plotBars, plotRasters, ...
                    plotMaps, summaryFig, plotLFP, showFigures);
                
                % Plot waveforms?
                if plotWFs
                    disp('Plotting waveforms...');
                    plotWaveforms([], figuresPath, fileName, Electrodes);
                end
                
                % Plot ISIs?
                if plotISI && ~isempty(Results)
                    disp('Plotting ISIs...');
                    plotISIs(figuresPath, fileName, Results);
                end
        end
    otherwise
end

