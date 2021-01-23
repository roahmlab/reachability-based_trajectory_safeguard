classdef rlSGDMSolver < rl.internal.optimizer.rlAbstractSGDMSolver
    % RLSOLVERSGDM   Stochastic gradient descent with momentum (SGDM) solver.
    % (wrapper of sgdmupdate).
    
    %   Copyright 2019 The MathWorks, Inc.
    
    methods
        function this = rlSGDMSolver(Options)
            % Constructor
            
            this = this@rl.internal.optimizer.rlAbstractSGDMSolver(Options);
        end
        
        function [this, LearnableParameters] = calculateUpdate(this,LearnableParameters,Gradients,GlobalLearnRate)
            [LearnableParameters, this.PreviousVelocities] = sgdmupdate(LearnableParameters, ...
                Gradients, this.PreviousVelocities, GlobalLearnRate, this.Momentum);
        end
    end
end