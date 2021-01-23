%rlValueRepresentation State value representation V(o).

% Copyright 2019 The MathWorks, Inc.

classdef rlValueRepresentation < rl.representation.rlAbstractValueRepresentation
    
    methods
        function this = rlValueRepresentation(Model, ObservationInfo, Options)
            % Constructor
            ActionInfo = {};
            this = this@rl.representation.rlAbstractValueRepresentation(Model, ObservationInfo, ActionInfo, Options);
            
            % validate model input dimensions
            % - check if the model has the same number of input channel
            %   specified by obs specs
            % - check if each output channel of the model has compatible
            %   size with obs specs dimension
            validateModelInputDimension(this)
            
            % validate model input dimensions: Vrep has single output
            ModelOutputSize = getSize(this.Model, 'output');
            NumOutput = prod(ModelOutputSize{1});
            if NumOutput ~= 1
                error(message('rl:agent:errStateValueRepNotSingleOutput'));
            end
            
            % default to use mse loss
            this = setLoss(this,rl.util.getDefaultValueRepLoss());
        end
    end
    
    methods (Access = protected)
        function [StateValue, State] = getValueImpl(this, Observation)
            % Compute V(o) - state value estimation of a given observation
            % Input:
            %   Observation: cell array for each observation channel.
            %       Each cell contains observation in obsDim x nBatch x nSequence
            % Output:
            %   StateValue: numeric array 1 x nBatch x nSequence
            %   State: State of previous prediction (e.g RNN hidden state)
            
            [StateValue, State, BatchSize, SequenceLength] = evaluate(this, Observation);
            % Output in format 1 x BatchSize x SequenceLength
            StateValue = reshape(StateValue{:},[1 BatchSize SequenceLength]);
        end
    end
    
end