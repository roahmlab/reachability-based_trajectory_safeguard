%rlStochasticActorRepresentation Stochastic actor representation pi(o).

% Copyright 2019 The MathWorks, Inc.

classdef rlStochasticActorRepresentation < rl.representation.rlAbstractActorRepresentation
    % REVISIT: support multiple action channels
    % REVISIT: support other/custom sampling strategy
    
    properties (Hidden,SetAccess = private)
        % TODO: make private and introduce method to compute pdf
        % Stochastic policy sampling strategy. Specify how to sample action
        % from probability parameters and the probability function
        SamplingStrategy
    end
    
    methods
        function this = rlStochasticActorRepresentation(Model, ObservationInfo, ActionInfo, Options)
            this = this@rl.representation.rlAbstractActorRepresentation(Model, ObservationInfo, ActionInfo, Options);
            
            % only support single channel outputs
            OutputSize = getSize(this,'output');
            
            if rl.util.isaSpecType(this.ActionInfo, 'continuous')
                % only support vector continuous action
                if any(arrayfun(@(x) iCheckNotVectorAction(x), this.ActionInfo))
                    error(message('rl:agent:errContinuousStochasticActorRepActionNotVector'))
                end
                
                % number of output must be 2 * number of actions
                NumContAction = sum(prod(this.ActionInfo.Dimension),'all');
                if prod(OutputSize{1}) ~= 2*NumContAction
                    error(message('rl:agent:errContinuousStochasticActorRepNumOutNeq2xNumAct'))
                end
                
                % set up sampling strategy
                this.SamplingStrategy = rl.representation.sampling.rlContinuousGaussianStrategy(this.ActionInfo);
                
                % set up model to gaussian mean and standard deviation
                this.Model = setModelOutputType(this.SamplingStrategy, this.Model);
            else
                % number of outputs must be equal to number of discrete actions
                NumDiscreteAction = prod(getNumberOfElements(this.ActionInfo));
                if prod(OutputSize{1}) ~= NumDiscreteAction
                    error(message('rl:agent:errDiscreteStochasticActorRepNumOutNeqNumAct'))
                end
                
                % set up sampling strategy
                this.SamplingStrategy = rl.representation.sampling.rlDiscreteSamplingStrategy(this.ActionInfo);
                
                % set up model to output probability
                this.Model = setModelOutputType(this.SamplingStrategy, this.Model);
            end
            
            % TODO: validate. Should sampling strategy validation be here?
            validateModelOutputCompatibility(this.SamplingStrategy, getSize(this.Model,'output'));
        end
        
        function Policy = evaluatePolicy(this, ActionProbabilityParam, Action)
            % Evaluate the policy density function. e.g. pi(x|mu,sigma)
            % ActionProbabilityParam (e.g. mean, standard deviation) at Action (x)
            %   Policy = ActionProbabilityParam (action probabilities) for 
            %   discrete action
            
            Policy = evaluate(this.SamplingStrategy, ActionProbabilityParam, Action);
        end
    end
    
    methods (Access = protected)
        function [Action, State, BatchSize, SequenceLength] = getActionImpl(this, Observation)
            % Return the action sampled from the current policy and
            % observation
            
            [ProbabilityParameter, State, BatchSize, SequenceLength] = evaluate(this, Observation);
            Action = sample(this.SamplingStrategy, ProbabilityParameter, BatchSize, SequenceLength);
        end
    end
end

function tf = iCheckNotVectorAction(ActionInfo)

tf = (numel(ActionInfo.Dimension) > 2) && all(ActionInfo.Dimension ~= 1);
end