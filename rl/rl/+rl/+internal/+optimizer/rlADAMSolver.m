classdef rlADAMSolver < rl.internal.optimizer.rlAbstractADAMSolver
    % RLSOLVERADAM   Adaptive Moment estimation (ADAM) solver (wrapper of 
    % adamupdate).
    
    %   Copyright 2019 The MathWorks, Inc.
    
    methods
        function this = rlADAMSolver(Options)
            % Constructor
            
            this = this@rl.internal.optimizer.rlAbstractADAMSolver(Options);
        end
        
        function [this, LearnableParameters] = calculateUpdate(this,LearnableParameters,Gradients,GlobalLearnRate)
            this.NumUpdates = this.NumUpdates + 1;
            [LearnableParameters, this.GradientMovingAverage, this.SquaredGradientMovingAverage] = ...
                adamupdate(LearnableParameters, Gradients, ...
                this.GradientMovingAverage, this.SquaredGradientMovingAverage, this.NumUpdates, ...
                GlobalLearnRate, this.GradientDecayFactor, this.SquaredGradientDecayFactor, this.Epsilon);
        end
    end
end
