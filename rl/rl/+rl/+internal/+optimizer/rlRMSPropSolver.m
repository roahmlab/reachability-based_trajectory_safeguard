classdef rlRMSPropSolver < rl.internal.optimizer.rlAbstractRMSPropSolver
    % RLSOLVERRMSPROP   Root Mean Square Propagation (RMSProp) solver 
    % (wrapper of rmspropupdate).
    
    %   Copyright 2019 The MathWorks, Inc.
    
    methods
        function this = rlRMSPropSolver(Options)
            % Constructor
            
            this = this@rl.internal.optimizer.rlAbstractRMSPropSolver(Options);
        end
        
        function [this, LearnableParameters] = calculateUpdate(this,LearnableParameters,Gradients,GlobalLearnRate)            
            [LearnableParameters, this.SquaredGradientMovingAverage] = rmspropupdate(...
                LearnableParameters,Gradients,this.SquaredGradientMovingAverage,...
                GlobalLearnRate,this.SquaredGradientDecayFactor,this.Epsilon);
        end
    end
end