classdef AbstractQDelegate
    %ABSTRACTQDELEGATE Interface for Q delegator. Q delegator performs
    % error checking, different getValue and getMaxQValue for single-output
    % and multi-output Q representations.
    
    % Copyright 2019 The MathWorks, Inc.
    
    properties
        % Determines the output size of getValue. 
        %  OutputSize = 1 for single output Q
        %  OutputSize = NumDiscreteAction for multi output Q
        OutputSize
    end
    
    methods
        function [QValue, State] = getValue(this, QRep, varargin)
            
            
            [QValue, State, BatchSize, SequenceLength] = getValueImpl(this, QRep, varargin{:});
            % output in format OutputSize x BatchSize x SequenceLength
            QValue = reshape(QValue{:},[this.OutputSize BatchSize SequenceLength]);
        end
        
        function [MaxQ, MaxIndex, State] = getMaxQValue(this, QRep, Observation)
            
            [MaxQ, MaxIndex, State, BatchSize, SequenceLength] = getMaxQValueImpl(this, QRep, Observation);
            % output in format 1 x BatchSize x SequenceLength
            MaxQ = reshape(MaxQ, 1, BatchSize, SequenceLength);
            MaxIndex = reshape(MaxIndex, 1, BatchSize, SequenceLength);
        end
    end
    
    methods (Access = protected)
        [QValue, State, BatchSize, SequenceLength] = getValueImpl(this, QRep, Observation);
        [MaxQ, MaxIndex, State] = getMaxQValueImpl(this, QRep, Observation)
    end
end

