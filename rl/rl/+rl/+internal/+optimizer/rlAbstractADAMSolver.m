classdef rlAbstractADAMSolver < rl.internal.optimizer.rlAbstractSolver
    % RLSOLVERADAM   Adaptive Moment estimation (ADAM) solver (wrapper of 
    % adamupdate) for framework 2.0 dlnetwork.
    
    %   Copyright 2019 The MathWorks, Inc.
    
    %   This utility will be replaced with DLT FW2.0 external feature. 
    %   Borrow from nnet.internal.cnn.solver.SolverADAM
    
    properties
        % GradientDecayFactor   Decay factor for moving average of gradients.
        %   A real scalar in [0,1) specifying the exponential decay rate
        %   for the gradient moving average. This parameter is denoted by
        %   the symbol Beta1 in the ADAM paper.
        GradientDecayFactor
        
        % SquaredGradientDecayFactor   Decay factor for moving average of squared gradients.
        %   A real scalar in [0,1) specifying the exponential decay rate
        %   for the squared gradient moving average. This parameter is
        %   denoted by the symbol Beta2 in the ADAM paper.
        SquaredGradientDecayFactor
        
        % Epsilon   Offset for the denominator in the ADAM update.
        %   A positive real scalar specifying the offset to use in the
        %   denominator for the ADAM update to prevent divide-by-zero
        %   problems.
        Epsilon
    end
    
    properties(Access = protected)
        % GradientMovingAverage   Moving average of gradients.
        %   A cell array of length NumLearnableParameters. Each element of
        %   the cell array contains the moving average of gradient for that
        %   learnable parameter.
        GradientMovingAverage
        
        % SquaredGradientMovingAverage   Moving average of squared gradients.
        %   A cell array of length NumLearnableParameters. Each element of
        %   the cell array contains the moving average of the squared
        %   gradient for that learnable parameter.
        SquaredGradientMovingAverage
        
        % NumUpdates   Number of updates so far.
        %   A non-negative integer indicating the number of update steps
        %   that have been computed so far.
        NumUpdates
    end
    
    methods
        function this = rlAbstractADAMSolver(Options)
            % Constructor
            
            this.GradientDecayFactor = Options.OptimizerParameters.GradientDecayFactor;
            this.SquaredGradientDecayFactor = Options.OptimizerParameters.SquaredGradientDecayFactor;
            this.Epsilon = Options.OptimizerParameters.Epsilon;
            this.NumUpdates = 0;
        end
    end
    
    %======================================================================
    % Save/Load
    %======================================================================
    methods
        function this = saveobj(this)
            % Save statistics as cpu array so users can load from non-gpu
            % machines
            if ~isempty(this.GradientMovingAverage)
                UseGPU = rl.internal.optimizer.rlAbstractSolver.useGPU(this.GradientMovingAverage);
                if UseGPU
                    this.GradientMovingAverage = dlupdate(@gather,this.GradientMovingAverage);
                    this.SquaredGradientMovingAverage = dlupdate(@gather,this.SquaredGradientMovingAverage);
                end
            end
        end
    end
end
