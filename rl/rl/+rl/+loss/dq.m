function Loss = dq(QPrediction,LossVariable)
    % Deep Q-Learning loss.
    %   QPrediction: output from model
    %   LossVariable: struct contains any info necessary to compute loss
    
    ObsDim = finddim(QPrediction,'B');
    NumObs = size(QPrediction, ObsDim);
    
    % Bellman eqn
    TargetQValues = getMaxQValue(LossVariable.TargetCritic, LossVariable.NextObservation);
    TargetQValues(~LossVariable.DoneIdx) = LossVariable.Reward(~LossVariable.DoneIdx) + ...
        LossVariable.Discount.*TargetQValues(~LossVariable.DoneIdx);
    
    % for terminal step, use the immediate reward (no more next state)
    TargetQValues(LossVariable.DoneIdx) = LossVariable.Reward(LossVariable.DoneIdx);
    if strcmpi(getQType(LossVariable.TargetCritic), 'multiOutput')
        ActionIdxMat = getElementIndicationMatrix(LossVariable.ActionInfo,LossVariable.Action,NumObs);
        
        % REVISIT: Abstract dataformat for robustness
        QPrediction = QPrediction(ActionIdxMat);
        Loss = mse(QPrediction, reshape(TargetQValues,size(QPrediction)),'DataFormat','BC');
    else
        Loss = mse(QPrediction, reshape(TargetQValues,size(QPrediction)));
    end
end