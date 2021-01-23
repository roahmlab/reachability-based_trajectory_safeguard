function Loss = policyGradientContinuous(MeanAndStd, LossVariable)
% Policy gradient loss with weighted entropy loss to encourages policy
% exploration.
% TotalLoss = PolicyGradientLoss + Weight*EntropyLoss
%   MeanAndVariance: mean and std batched together - dlarray output from model
%   LossVariable: struct contains Advantage, previous Action, 
%                 Sampling strategy and EntropyLossWeight

% Copyright 2019 The MathWorks Inc.

NumExperience = numel(LossVariable.Advantage);

% policy gradient loss
Policy = evaluate(LossVariable.SamplingStrategy, MeanAndStd, LossVariable.Action);
LogPolicy = log(rl.internal.dataTransformation.boundAwayFromZero(Policy));
PolicyGradientLoss = -sum(LogPolicy .* LossVariable.Advantage, 'all')/NumExperience;

% entropy loss
EntropyLoss = rl.loss.policyEntropyContinuous(MeanAndStd, ...
    LossVariable.EntropyLossWeight,NumExperience);

% total loss 
Loss = PolicyGradientLoss + EntropyLoss;
end