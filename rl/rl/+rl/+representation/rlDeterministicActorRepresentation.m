%rlDeterministicActorRepresentation Deterministic actor representation pi(o).

% Copyright 2019 The MathWorks, Inc.

classdef rlDeterministicActorRepresentation < rl.representation.rlAbstractActorRepresentation
    methods
        function this = rlDeterministicActorRepresentation(Model, ObservationInfo, ActionInfo, Options)
            % Constructor
            
            this = this@rl.representation.rlAbstractActorRepresentation(Model, ObservationInfo, ActionInfo, Options);
            
            % deterministic actor not support discrete action
            if any(arrayfun(@(x) isa(x,'rl.util.rlFiniteSetSpec'), ActionInfo))
                error(message('rl:agent:errDetermActorNotSupportDiscreteAction'));
            end
            
            % validate if model size is compatible with action info
            validateActionOutputDimension(this)
        end
    end
    
    methods (Access = protected)
        function [Action, State, BatchSize, SequenceLength] = getActionImpl(this, Observation)
            % Return the deterministic action from the current policy and
            % observation
            
            [Action, State, BatchSize, SequenceLength] = evaluate(this, Observation);
        end
    end
    
    methods (Access = private)
        function validateActionOutputDimension(this)
            % Check if model size is compatible with action info
            
            ModelOutputSize = reshape(getSize(this.Model,'output'),1,[]);
            
            if ~all(cellfun(@(x,y) rl.util.isArrayEqual(x,y), ModelOutputSize, this.ActionDimension))
                error(message('rl:agent:errIncompatibleDetermActionDim'));
            end
        end
    end
end