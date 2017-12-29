classdef visstimconfig
    %VISSTIMCONFIG Configuration for VisStim
    %   Detailed explanation goes here
    
    properties(Constant)
        ParallelDCh=1;
        PhotodiodeDCh=2;
        StartDCh=3;
        StopDCh=4
        MarkSearchRadius=0.050; % s
        NMarkPerCond=2; % number of mark per condition
        TimerDriftError=[0;0.000036];
    end
    
    methods
    end
end

