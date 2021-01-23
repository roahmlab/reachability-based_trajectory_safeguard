classdef rlLayerADAMSolver < rl.internal.optimizer.rlAbstractADAMSolver & ...
                             rl.internal.optimizer.rlAbstractLayerSolver
    % RLSOLVERADAM   Adaptive Moment estimation (ADAM) solver for internal
    % layer API. This is a value object implementation of 
    % nnet.internal.cnn.solver.SolverADAM
    
    %   Copyright 2019 The MathWorks, Inc.
    
    %   This utility will be removed once RLT transit to dlnetwork
    
    methods
        function this = rlLayerADAMSolver(Options, Model)
            % Constructor
            
            this = this@rl.internal.optimizer.rlAbstractADAMSolver(Options);
            this = this@rl.internal.optimizer.rlAbstractLayerSolver(Model);
            this = initializeState(this);
            % cast to single to match nnet.internal.cnn.solver.SolverADAM
            this.GradientDecayFactor = single(this.GradientDecayFactor);
            this.SquaredGradientDecayFactor = single(this.SquaredGradientDecayFactor);
            this.Epsilon = single(this.Epsilon);
        end
        
        function [this, LearnableParameters] = calculateUpdate(this,LearnableParameters,Gradients,GlobalLearnRate)
            
            LocalLearnRates = this.LocalLearnRates;
            NumLearnableParameters = this.NumLearnableParameters;
            
            Beta1 = this.GradientDecayFactor;
            Beta2 = this.SquaredGradientDecayFactor;
            Epsilon = this.Epsilon;
            
            this.NumUpdates = this.NumUpdates + 1;
            LearnRateShrinkFactor = iCalculateLearnRateShrinkFactor(Beta1,Beta2,this.NumUpdates);
            
            for i = 1:NumLearnableParameters
                % No update needed for parameters that are not learning
                if any(LocalLearnRates{i},'all') && ~isempty(Gradients{i})
                    EffectiveLearningRate = LearnRateShrinkFactor.*GlobalLearnRate.*LocalLearnRates{i};
                    [Gradients{i}, this.GradientMovingAverage{i}, this.SquaredGradientMovingAverage{i}] = nnet.internal.cnn.solver.adamstep(...
                        Gradients{i}, this.GradientMovingAverage{i}, this.SquaredGradientMovingAverage{i}, EffectiveLearningRate, Beta1, Beta2, Epsilon);
                    
                    % add gradients to get updated LearnableParameters
                    % similar to internal Layer API updateLearnableParameters()
                    LearnableParameters{i} = LearnableParameters{i} + Gradients{i};
                end
            end
        end
    end
    
    methods (Access = protected)
        function this = initializeState(this)
            % initializeState(this) sets the state of the solver to its
            % initial state.
            
            this.GradientMovingAverage = iInitializeMovingAverage(this.NumLearnableParameters);
            this.SquaredGradientMovingAverage = iInitializeMovingAverage(this.NumLearnableParameters);
        end
    end
end

function MovingAverage = iInitializeMovingAverage(NumLearnableParameters)
MovingAverage = cell(1,NumLearnableParameters);
for i = 1:NumLearnableParameters
    MovingAverage{i} = single(zeros(1));
end
end

function LearnRateShrinkFactor = iCalculateLearnRateShrinkFactor(Beta1,Beta2,t)
LearnRateShrinkFactor = sqrt(1-Beta2^t)/(1-Beta1^t);
end