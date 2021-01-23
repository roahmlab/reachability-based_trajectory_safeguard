classdef rlAbstractLayerSolver
    % Solver   Abstract class for internal layer API solvers.
    
    %   Copyright 2019 The MathWorks, Inc.
    
    % Will be removed when RLT transits to dlnetwork
    
    properties (Access = protected)
        % NumLearnableParameters   Number of learnable parameters.
        %   A non-negative integer specifying the number of learnable
        %   parameters in the optimization problem.
        NumLearnableParameters
        
        % LocalLearnRates   Local learning rates for learnable parameters.
        %   A cell array of length NumLearnableParameters where each
        %   element of the cell array specifies the local learning rate for
        %   that learnable parameter. Every learnable parameter can define
        %   its own local learning rate. The effective learning rate for
        %   each learnable parameter is calculated by taking the product of
        %   the current global learning rate and its local learning rate.
        LocalLearnRates
    end
    
    methods
        function this = rlAbstractLayerSolver(LayerModel)
            % this = Solver(learnableParameters,precision) creates a Solver
            % object for optimizing parameters in learnableParameters using
            % floating point precision specified in precision. The class of
            % inputs is as follows:
            %
            %    learnableParameters - an array of objects of type nnet.internal.cnn.layer.learnable.LearnableParameter
            %    precision           - an object of type nnet.internal.cnn.util.Precision
            
            % get internal layer DAG object
            InternalNetwork = getInternalModel(LayerModel);
            this.NumLearnableParameters = numel(InternalNetwork.LearnableParameters);
            this.LocalLearnRates = iExtractLocalLearnRates(InternalNetwork.LearnableParameters);
        end
    end
    
    methods(Abstract, Access = protected)
        initializeState(this);
    end
end

function LocalLearnRates = iExtractLocalLearnRates(LearnableParameters)
NumLearnableParameters = numel(LearnableParameters);
LocalLearnRates = cell(1,NumLearnableParameters);
for i = 1:NumLearnableParameters
    LocalLearnRates{i} = single(LearnableParameters(i).LearnRateFactor);
end
end