%rlAbstractQRepresentation Defines interface for state-action value 
% representation Q(o) and Q(o,a).

% Copyright 2019 The MathWorks, Inc.

classdef rlAbstractQRepresentation < rl.representation.rlAbstractValueRepresentation
    
    methods
        function this = rlAbstractQRepresentation(Model, varargin)
            % Constructor
            
            this = this@rl.representation.rlAbstractValueRepresentation(Model, varargin{:});
            
            % Q rep requires action info
            if isempty(this.ActionInfo)
                error(message('rl:agent:errQValueRepActionInfoRequired'));
            end
        end
        
        function [MaxQ, MaxActionIndex, State] = getMaxQValue(this, Observation)
            % GETMAXQVALUE return the estimated maximum state-action value 
            % of the discrete action Q representation given the observation.
            %
            %   [MAXQ, MAXACTIONINDEX, STATE] = GETMAXQVALUE(QVALUEREP, OBS)
            %   returns the estimated maximum state-action value MAXQ of
            %   the Q value representation QVALUEREP at observation OBS.
            %
            %       OBS is a cell array with as many elements as the number
            %       of representation input channels. Each cell array element 
            %       must follow the dimension convention 
            %       SingleDataDimension-by-BatchSize-by-SequenceLength
            %
            %       MAXQ is a 1-by-BatchSize-by-SequenceLength matrix.
            %       MAXACTIONINDEX is a 1-by-BatchSize-by-SequenceLength 
            %       matrix of action index that results in MAXQ.
            
            % only support discrete action
            if isa(this.ActionInfo, 'rl.util.rlNumericSpec')
                error(message('rl:agent:errGetMaxQContinuousAction'));
            end
            [MaxQ, MaxActionIndex, State] = getMaxQValueImpl(this, Observation);
        end
    end
    
    methods (Abstract, Access = protected)
        [MaxQ, MaxIndex, State] = getMaxQValueImpl(this, Observation)
    end
    
end