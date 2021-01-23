function Loss = drq(QPrediction,LossVariable)
    % Deep Recurrent Q-Learning loss.
    %   QPrediction: output from model
    %   LossVariable: struct contains the following info necessary to compute loss
    %   TargetCritic: [1×1 rl.representation.rlQValueRepresentation]
    %       Action: {[Dimension x Batch size × sequenceLength]} single
    %       Reward: [1 × Batch size × sequenceLength] single
    %       NextObservation: {[Dimension × Batch size × sequenceLength]}  single
    %       DoneIdx: [1 × Batch size × sequenceLength] logical
    %       Discount: 1 × 1
    %       ActionInfo: [1×1 rl.util.rlFiniteSetSpec]
    %       MaskIdx: [1 × Batch size × sequenceLength]  single
    %       MasIkdx is zero if the corresponding input is the padded input.
    %
    
    BatchDimIdx =finddim(QPrediction,'B');
    TrajectoryDimIdx = finddim(QPrediction,'T');
    MiniBatchSize = size(QPrediction, BatchDimIdx);
    SequenceLength = size(QPrediction, TrajectoryDimIdx);
    % Save the data format of QPrediction. This is used when loss is
    % computed
    DataFormatQprediction = dims(QPrediction);
    % Compute critic target
    LossVariable.TargetCritic = resetState(LossVariable.TargetCritic);
    
    % Bellman eqn
    TargetQValues = getMaxQValue(LossVariable.TargetCritic, LossVariable.NextObservation);
    TargetQValues(~LossVariable.DoneIdx) = LossVariable.Reward(~LossVariable.DoneIdx) + ...
        LossVariable.Discount.*TargetQValues(~LossVariable.DoneIdx);
    
    % for terminal step, use the immediate reward (no more next state)
    TargetQValues(LossVariable.DoneIdx) = LossVariable.Reward(LossVariable.DoneIdx);
    
    if strcmpi(getQType(LossVariable.TargetCritic), 'multiOutput')
        
        % LossVariable.Action contains raw actions values. getElementIndicationMatrix first
        % converts them to integers (1,2,...,|A|), then create one-hot
        % encoding vectors.
        ActionIdxMat = getElementIndicationMatrix(LossVariable.ActionInfo,LossVariable.Action,MiniBatchSize*SequenceLength);
        ActionIdxMat = reshape(ActionIdxMat,[],MiniBatchSize,SequenceLength);
        
        % Compute Q values of selected actions.
        QPrediction = QPrediction(ActionIdxMat);
        
        % To ignore the effect of the padded inputs (MasIkdx of the padded samples are false),
        % we want to make the corresponding gradients to be zero. In order to do it, we use the
        % dummy Q values.
        
        % 1. QPrediction is dlarray. We convert the type of TargetQvalues
        % to be dlarray as well.
        TargetQValues = dlarray(TargetQValues);
        
        % 2. We use dummy Q values for the padded samples so that the
        % targets and q values are the same, hence there are no gradients
        % for the samples.
        
        % When it uses GPU for critic, the type of QPrediction is dlarray gpuArray.
        % We first move TargetQValiues to gpu. Then, it copies dummy Q
        % Values.
        if strcmp(LossVariable.TargetCritic.Options.UseDevice,"gpu")
            % Move TargetQValues to gpu memory
            TargetQValues = gpuArray(TargetQValues);
            % copy dummy Q values
            TargetQValues(~LossVariable.MaskIdx) = QPrediction(~LossVariable.MaskIdx); 
        else
            % copy dummy Q values
            TargetQValues(~LossVariable.MaskIdx) = QPrediction(~LossVariable.MaskIdx); 
        end
        
        % Compute loss        
        QPrediction = reshape(QPrediction, size(TargetQValues));
        Loss = mse(QPrediction, TargetQValues, 'DataFormat',DataFormatQprediction); %'CBT'       
        
    else
        % REVISIT after it supports multiple inputs for LSTM. It currently
        % does not use MaskIdx here.
        error((message('rl:agent:errStatefulSingleOutQDiscreteAction')));
        Loss = mse(QPrediction, reshape(TargetQValues,size(QPrediction)));
    end
end