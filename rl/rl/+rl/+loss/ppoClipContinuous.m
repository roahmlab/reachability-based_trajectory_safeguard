function Loss = ppoClipContinuous(MeanAndStd, LossVariable)
% NOTE: will replace rl.loss.ppoClipped
% Clipped PPO with entropy loss function function for continuous action space
%   MeanAndStd: dlarray of current policy action probabilities (model output)
%   LossVariable: struct contains
%       - SamplingStrategy
%       - Action: previous action
%       - OldPolicy: old action policy piOld(at|st)
%       - Advantage
%       - ClipFactor: scalar > 0
%       - EntropyLossWeight: scalar where 0 <= EntropyLossWeight <= 1

% Copyright 2019 The MathWorks Inc.

% Extract information from input
Advantage = LossVariable.Advantage;
OldPolicy = LossVariable.OldPolicy;
NumExperience = numel(Advantage);

% compute pi(at|st)
Policy = evaluate(LossVariable.SamplingStrategy, MeanAndStd, LossVariable.Action);
% rt = pi(at|st)/piOld(at|st), avoid division by zero
Ratio = Policy ./ rl.internal.dataTransformation.boundAwayFromZero(OldPolicy);

% obj = rt * At
Advantage = reshape(Advantage, 1, NumExperience);
Objective = Ratio .* Advantage;
ObjectiveClip = max(min(Ratio, 1 + LossVariable.ClipFactor), 1 - LossVariable.ClipFactor) .* Advantage;

% clipped surrogate loss
SurrogateLoss = -sum(min(Objective, ObjectiveClip),'all')/NumExperience;

% entropy loss
EntropyLoss = rl.loss.policyEntropyContinuous(MeanAndStd, ...
    LossVariable.EntropyLossWeight,NumExperience);

% total loss
Loss = SurrogateLoss + EntropyLoss;
end