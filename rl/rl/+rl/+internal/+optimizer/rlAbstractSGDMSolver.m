classdef rlAbstractSGDMSolver < rl.internal.optimizer.rlAbstractSolver
    % RLSOLVERSGDM   Stochastic gradient descent with momentum (SGDM) solver.
    % (wrapper of sgdmupdate) for framework 2.0 dlnetwork.
    
    %   Copyright 2019 The MathWorks, Inc.
    
    %   This utility will be replaced with DLT FW2.0 external feature. 
    %   Borrow from nnet.internal.cnn.solver.SolverSGDM
    
    properties
        % Momentum   Momentum for SGDM update.
        %   A real scalar in [0,1] specifying the coefficient of the
        %   momentum term in SGDM update.
        Momentum
    end
    
    properties (Access = protected)
        % PreviousVelocities   Previous update steps in SGDM.
        %   A cell array of length NumLearnableParameters. Each element of
        %   the cell array contains the update step computed in the
        %   previous iteration of SGDM for that learnable parameter.
        PreviousVelocities
    end
    
    methods
        function this = rlAbstractSGDMSolver(Options)
            % Constructor
            
            this.Momentum = Options.OptimizerParameters.Momentum;
        end
    end
    
    %======================================================================
    % Save/Load
    %======================================================================
    methods
        function this = saveobj(this)
            % Save statistics as cpu array so users can load from non-gpu
            % machines.
            if ~isempty(this.PreviousVelocities)
                UseGPU = rl.internal.optimizer.rlAbstractSolver.useGPU(this.PreviousVelocities);
                if UseGPU
                    this.PreviousVelocities = dlupdate(@gather,this.PreviousVelocities);
                end
            end
        end
    end
end