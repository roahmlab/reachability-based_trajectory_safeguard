%rlAbstractValueRepresentation Defines the interface for value functions.
% Value functions can be V(o), Q(o,a) and Q(o)

% Copyright 2019 The MathWorks, Inc.

classdef rlAbstractValueRepresentation < rl.representation.rlAbstractRepresentation
    
    methods
        function this = rlAbstractValueRepresentation(Model, ObservationInfo, ActionInfo, Options)
            % Constructor
            this = this@rl.representation.rlAbstractRepresentation(Model, ObservationInfo, ActionInfo, Options);
            
            % V, Q representation must have single channel output
            ModelOutputSize = getSize(this.Model, 'output');
            if numel(ModelOutputSize) > 1
                % REVISIT MSG replace 'channel' to be more user-friendly
                error(message('rl:agent:errValueRepNotSingleChannel'));
            end
        end
        
        function [Value, State] = getValue(this, varargin)
            % GETVALUE return the estimated state value or state-action 
            % value of the representation given the observation or 
            % observation-action pair.
            %
            %   [VALUE, STATE] = GETVALUE(VALUEREP, OBS) returns
            %   the estimated state value VALUE and current state STATE of 
            %   the state value representation VALUEREP at observation OBS.
            %
            %       OBS is a cell array with as many elements as the number
            %       of representation input channels. Each cell array element 
            %       must follow the dimension convention 
            %       SingleDataDimension-by-BatchSize-by-SequenceLength
            %
            %       VALUE is a 1-by-BatchSize-by-SequenceLength matrix.
            %
            %   [VALUE, STATE] = GETVALUE(QVALUEREP, OBS) returns the
            %   estimated state-action value VALUE and current state STATE 
            %   of the Q value representation QVALUEREP at observation OBS.
            %   QVALUEREP must have as many outputs as the number of
            %   discrete actions.
            %
            %       VALUE is a NumAction-by-BatchSize-by-SequenceLength 
            %       matrix.
            %
            %   [VALUE, STATE] = GETVALUE(QVALUEREP, OBS, ACT) returns the
            %   estimated state-action value VALUE and current state STATE 
            %   of the Q value representation QVALUEREP at observation OBS
            %   and action ACT. QVALUEREP must have single output.
            %
            %       OBS and ACT are cell arrays with as many elements as 
            %       the number of observation and action channels. Each 
            %       cell array element must follow the dimension convention 
            %       SingleDataDimension-by-BatchSize-by-SequenceLength
            %
            %       VALUE is a 1-by-BatchSize-by-SequenceLength 
            %       matrix.
            
            [Value, State] = getValueImpl(this, varargin{:});
            % return numeric value (always single channel)
            if iscell(Value)
                Value = Value{:};
            end
        end
    end
    
    methods (Abstract, Access = protected)
        [Value, State] = getValueImpl(this, varargin)
    end
    
end