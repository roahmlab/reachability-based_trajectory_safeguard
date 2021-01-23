classdef QMultiOutputDelegate < rl.representation.qdelegate.AbstractQDelegate
    % Multi output Q delegator performs error checking, getValue and
    % getMaxQValue for multi-output Q representations.
    
    % Copyright 2019 The MathWorks, Inc.

    methods
        function this = QMultiOutputDelegate(QRep)
            
            % check if ActionInfo is discrete
            iValidateDiscreteAction(QRep)
            
            % check if number of output is equal to number of discrete action
            NumDiscreteAction = iValidateNumOutputNumAction(QRep);
            this.OutputSize = NumDiscreteAction;
        end
    end
    
    methods (Access = protected)
        function [QValue, State, BatchSize, SequenceLength] = getValueImpl(~, QRep, Observation)
            [QValue, State, BatchSize, SequenceLength] = evaluate(QRep, Observation);
        end

        function [MaxQ, MaxIndex, State, BatchSize, SequenceLength] = getMaxQValueImpl(~, QRep, Observation)
            
            [QValue, State, BatchSize, SequenceLength] = evaluate(QRep, Observation);
            [MaxQ, MaxIndex] = max(QValue{:});
        end
    end
end

function iValidateDiscreteAction(QRep)
% check if ActionInfo is discrete
if rl.util.isaSpecType(QRep.ActionInfo, 'continuous')
    error(message('rl:agent:errQContinuousActionMultiOutput'));
end
end

function NumDiscreteAction = iValidateNumOutputNumAction(Qrep)
% check if number of output is equal to number of discrete action

NumDiscreteAction = prod(getNumberOfElements(Qrep.ActionInfo),'all');
ModelOutputSize = getSize(Qrep, 'output');
NumOutput = prod(ModelOutputSize{1});
if NumDiscreteAction ~= NumOutput
    error(message('rl:agent:errQMultiNumOutputNotEqualNumDiscreteAction'));
end

end