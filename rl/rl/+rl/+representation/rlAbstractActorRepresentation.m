%rlAbstractActorRepresentation Defines the interface for actor (policy).

% Copyright 2019 The MathWorks, Inc.

classdef rlAbstractActorRepresentation < rl.representation.rlAbstractRepresentation
    
    properties (Dependent, SetAccess = private)
        ActionInfo
    end
    
    properties (Access = protected)
        % Action data dimension to enforce getAction data convention
        ActionDimension
    end
    
    methods
        function this = rlAbstractActorRepresentation(Model, varargin)
            % Constructor
            this = this@rl.representation.rlAbstractRepresentation(Model, varargin{:});
            
            if isa(this.Model,'rl.representation.model.rlTableModel')
                error(message('rl:agent:errTableActorNotSupport'))
            end
            
            % actor rep requires action info
            if isempty(this.ActionInfo)
                error(message('rl:agent:errActorRepActionInfoRequired'));
            end
            
            % action data dimension to enforce getAction data convention
            this.ActionDimension = {this.ActionInfo.Dimension};
            
            % validate model input dimensions
            % - check if the model has the same number of input channel
            %   specified by obs specs
            % - check if each output channel of the model has compatible
            %   size with obs specs dimension
            validateModelInputDimension(this)
        end
        
        function [Action, State] = getAction(this, Observation)
            % GETACTION return the action derived from the actor (policy) 
            % given the observation.
            %
            %   [ACTION, STATE] = GETACTION(ACTORREP, OBS) returns
            %   the action ACTION and current state STATE derived from the 
            %   actor representation ACTORREP at observation OBS.
            %
            %       OBS is a cell array with as many elements as the number
            %       of representation input channels. Each cell array element 
            %       must follow the dimension convention 
            %       ObservationDimension-by-BatchSize-by-SequenceLength
            %
            %       ACTION is a cell array with as many elements as the
            %       number of action channels. Each cell array element 
            %       follows the dimension convention 
            %       ActionDimension-by-BatchSize-by-SequenceLength
            %       If the ACTORREP has single action channel, ACTION is
            %       numeric.
            
            [Action, State, BatchSize, SequenceLength] = getActionImpl(this, Observation);
            Action = reshapeActionOutput(this, Action, BatchSize, SequenceLength);
        end
        
        function ActionInfo = get.ActionInfo(this)
            ActionInfo = this.ActionInfo_;
        end
    end
    
    methods (Abstract, Access = protected)
        [Value, State] = getActionImpl(this, Observation)
    end
    
    methods (Access = protected)
        function Action = reshapeActionOutput(this, Action, BatchSize, SequenceLength)
            % Return cell array. Each element has form [ActionInfoDim BatchSize SequenceLength]
            
            Action = cellfun(@(x,y) reshape(x,[y BatchSize SequenceLength]), Action, this.ActionDimension, 'UniformOutput', false);
        end
    end
    
end