classdef QSingleOutputDelegate < rl.representation.qdelegate.AbstractQDelegate
    % Single output Q delegator performs error checking, getValue and
    % getMaxQValue for single-output Q representations.
    
    % Copyright 2019 The MathWorks, Inc.
    
    properties (Access = private)
        % Flag determines whether model uses dlarray
        UseDLArray
    end
    
    methods
        function this = QSingleOutputDelegate(QRep)
            
            % Check if RNN detected for discrete action. Only support multi output discrete Q RNN
            iValidateNotStatefulRepDiscreteAction(QRep)
            
            this.UseDLArray = isa(QRep.getRLModel,'rl.representation.model.rlBasisFunctionModel');
            
            this.OutputSize = 1;
        end
    end
    
    methods (Access = protected)
        function [QValue, State, BatchSize, SequenceLength] = getValueImpl(~, QRep, Observation, Action)
            
            % NOTE: include due to dlarray cell concat inconsistent g2153798
            if ~iscell(Observation)
                Observation = {Observation};
            end
            if ~iscell(Action)
                Action = {Action};
            end
            [QValue, State, BatchSize, SequenceLength] = evaluate(QRep, [Observation, Action]);
        end
        
        function [MaxQ, MaxIndex, State, BatchSize, SequenceLength] = getMaxQValueImpl(this, QRep, Observation)
            % REVISIT: support multiple action channels
            
            ActionInfo = QRep.ActionInfo;
            
            % Not support RNN for single output Q
            [BatchSize,SequenceLength] = inferDataDimension(QRep.ObservationInfo(1), Observation{1});
            if SequenceLength > 1
                error(message('rl:agent:errStatefulSingleOutQDiscreteAction'));
            end
            
            NumAction = getNumberOfElements(ActionInfo);
            % REVISIT: if we decide output of evaluate is gpuArray
            QValue = zeros(NumAction,BatchSize,'single');
            if this.UseDLArray
                QValue = dlarray(QValue);
            end
            
            for ct = 1:NumAction
                CurrentAction = getElementValue(ActionInfo,ct);
                Q = evaluate(QRep,[Observation,{repmat(CurrentAction,[ones(1,numel(ActionInfo.Dimension)) BatchSize])}]);
                QValue(ct,:) = Q{:};
            end
            
            [MaxQ, MaxIndex] = max(QValue,[],1);
            State = [];
            SequenceLength = 1;
        end
    end
end

function iValidateNotStatefulRepDiscreteAction(QRep)
% Check if all actions are inputs to the representation

if hasState(QRep) && isa(QRep.ActionInfo, 'rl.util.rlFiniteSetSpec')
    error(message('rl:agent:errStatefulSingleOutQDiscreteAction'));
end
end