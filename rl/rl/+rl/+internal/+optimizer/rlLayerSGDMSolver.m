classdef rlLayerSGDMSolver < rl.internal.optimizer.rlAbstractSGDMSolver & ...
                             rl.internal.optimizer.rlAbstractLayerSolver
    % SolverSGDM   Stochastic gradient descent with momentum (SGDM) solver.
    % This is a value object implementation of 
    % nnet.internal.cnn.solver.SolverSGDM
    
    %   Copyright 2017-2019 The MathWorks, Inc.
    
    %   This utility will be removed once RLT transit to dlnetwork
    
    methods
        function this = rlLayerSGDMSolver(Options, Model)
            % Constructor
            
            this = this@rl.internal.optimizer.rlAbstractSGDMSolver(Options);
            this = this@rl.internal.optimizer.rlAbstractLayerSolver(Model);
            this = initializeState(this);
            % cast to single to match nnet.internal.cnn.solver.SolverSGDM
            this.Momentum = single(this.Momentum);
        end
        
        function [this, LearnableParameters] = calculateUpdate(this,LearnableParameters,Gradients,GlobalLearnRate)
            
            LocalLearnRates = this.LocalLearnRates;
            NumLearnableParameters = this.NumLearnableParameters;
            Momentum = this.Momentum;
            
            for i = 1:NumLearnableParameters
                % No update needed for parameters that are not learning
                if any(LocalLearnRates{i},'all') && ~isempty(Gradients{i})
                    EffectiveLearningRate = GlobalLearnRate.*LocalLearnRates{i};
                    [Gradients{i}, this.PreviousVelocities{i}] = nnet.internal.cnn.solver.sgdmstep(...
                         Gradients{i}, this.PreviousVelocities{i}, EffectiveLearningRate, Momentum);
                     
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
            
            this.PreviousVelocities = iInitializeVelocities(this.NumLearnableParameters);
        end
    end
end

function Velocities = iInitializeVelocities(NumLearnableParameters)
Velocities = cell(1,NumLearnableParameters);
for i = 1:NumLearnableParameters
    Velocities{i} = single(zeros(1));
end
end