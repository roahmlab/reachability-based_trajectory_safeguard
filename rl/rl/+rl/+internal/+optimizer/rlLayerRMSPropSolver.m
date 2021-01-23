classdef rlLayerRMSPropSolver < rl.internal.optimizer.rlAbstractRMSPropSolver & ...
                                rl.internal.optimizer.rlAbstractLayerSolver
    % rlLayerRMSPropSolver   Root Mean Square Propagation (RMSProp) solver.
    % This is a value object implementation of nnet.internal.cnn.solver.SolverRMSProp
    %   Copyright 2019 The MathWorks, Inc.
    
    %   This utility will be removed once RLT transit to dlnetwork
    
    methods
        function this = rlLayerRMSPropSolver(Options, Model)
            % Constructor
            
            this = this@rl.internal.optimizer.rlAbstractRMSPropSolver(Options);
            this = this@rl.internal.optimizer.rlAbstractLayerSolver(Model);
            this = initializeState(this);
            % cast to single to match nnet.internal.cnn.solver.SolverRMSProp
            this.SquaredGradientDecayFactor = single(this.SquaredGradientDecayFactor);
            this.Epsilon = single(this.Epsilon);
        end
        
        function [this, LearnableParameters] = calculateUpdate(this,LearnableParameters,Gradients,GlobalLearnRate)
            
            LocalLearnRates = this.LocalLearnRates;
            NumLearnableParameters = this.NumLearnableParameters;
            
            Rho = this.SquaredGradientDecayFactor;
            Epsilon = this.Epsilon;
            
            for i = 1:NumLearnableParameters
                % No update needed for parameters that are not learning
                if any(LocalLearnRates{i},'all') && ~isempty(Gradients{i})
                    EffectiveLearningRate = GlobalLearnRate.*LocalLearnRates{i};
                    [Gradients{i}, this.SquaredGradientMovingAverage{i}] = nnet.internal.cnn.solver.rmspropstep(...
                        Gradients{i}, this.SquaredGradientMovingAverage{i}, EffectiveLearningRate, Rho, Epsilon);
                    
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
            
            this.SquaredGradientMovingAverage = iInitializeMovingAverage(this.NumLearnableParameters);
        end
    end
end

function MovingAverage = iInitializeMovingAverage(numLearnableParameters)
MovingAverage = cell(1,numLearnableParameters);
for i = 1:numLearnableParameters
    MovingAverage{i} = single(zeros(1));
end
end