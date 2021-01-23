classdef rlAbstractRMSPropSolver < rl.internal.optimizer.rlAbstractSolver
    % RLSOLVERRMSPROP   Root Mean Square Propagation (RMSProp) solver 
    % (wrapper of rmspropupdate) for framework 2.0 dlnetwork.
    
    %   Copyright 2019 The MathWorks, Inc.
    
    %   This utility will be replaced with DLT FW2.0 external feature. 
    %   Borrow from nnet.internal.cnn.solver.SolverRMSProp
    
    properties
        % SquaredGradientDecayFactor   Decay factor for moving average of squared gradients.
        %   A real scalar in [0,1) specifying the exponential decay rate
        %   for the squared gradient moving average.
        SquaredGradientDecayFactor
        
        % Epsilon   Offset for the denominator in the RMSProp update.
        %   A positive real scalar specifying the offset to use in the
        %   denominator for the RMSProp update to prevent divide-by-zero
        %   problems.
        Epsilon
    end
    
    properties(Access = protected)
        % SquaredGradientMovingAverage   Moving average of squared gradients.
        %   A cell array of length NumLearnableParameters. Each element of
        %   the cell array contains the moving average of the squared
        %   gradient for that learnable parameter.
        SquaredGradientMovingAverage
    end
    
    methods
        function this = rlAbstractRMSPropSolver(Options)
            % Constructor
            
            this.SquaredGradientDecayFactor = Options.OptimizerParameters.SquaredGradientDecayFactor;
            this.Epsilon = Options.OptimizerParameters.Epsilon;
        end
    end
    
    %======================================================================
    % Save/Load
    %======================================================================
    methods
        function this = saveobj(this)
            % Save statistics as cpu array so users can load from non-gpu
            % machines
            if ~isempty(this.SquaredGradientMovingAverage)
                UseGPU = rl.internal.optimizer.rlAbstractSolver.useGPU(this.SquaredGradientMovingAverage);
                if UseGPU
                    this.SquaredGradientMovingAverage = dlupdate(@gather,this.SquaredGradientMovingAverage);
                end
            end
        end
    end
end