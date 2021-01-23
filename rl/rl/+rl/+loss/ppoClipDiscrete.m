function Loss = ppoClipDiscrete(Policy, LossVariable)
% NOTE: will replace rl.loss.ppoClipped
% Clipped PPO with entropy loss function function for discrete action space
%   ModelOutput: dlarray of current policy action probabilities
%   LossVariable: struct contains
%       - Action: 
%           continuous action space: previous action
%           discrete action space: logical matrix of action index
%       - OldPolicy: old action policy piOld(at|st)
%       - Advantage
%       - ClipFactor: scalar > 0
%       - EntropyLossWeight: scalar where 0 <= EntropyLossWeight <= 1

% Copyright 2019 The MathWorks Inc.

% Extract information from input
Advantage = LossVariable.Advantage;
OldPolicy = LossVariable.OldPolicy;
ActionIdx = LossVariable.Action;
NumExperience = numel(Advantage);

% rt = pi(at|st)/piOld(at|st), avoid division by zero
Ratio = Policy(ActionIdx) ./ rl.internal.dataTransformation.boundAwayFromZero(OldPolicy(ActionIdx));
Ratio = reshape(Ratio,size(Advantage));

% obj = rt * At
Objective = Ratio .* Advantage;
ObjectiveClip = max(min(Ratio, 1 + LossVariable.ClipFactor), 1 - LossVariable.ClipFactor) .* Advantage;

% clipped surrogate loss
SurrogateLoss = -sum(min(Objective, ObjectiveClip),'all')/NumExperience;

% entropy loss
EntropyLoss = rl.loss.policyEntropyDiscrete(Policy, ...
    LossVariable.EntropyLossWeight,NumExperience);

% total loss
Loss = SurrogateLoss + EntropyLoss;
end