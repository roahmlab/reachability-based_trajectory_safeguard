function Loss = policyGradientDiscrete(Policy, LossVariable)
% Policy gradient loss with weighted entropy loss to encourages policy
% exploration.
% TotalLoss = PolicyGradientLoss + Weight*EntropyLoss
%   Policy: pi(a|o) - dlarray output from model
%   LossVariable: struct contains Advantage, previous Action, 
%                 Sampling strategy and EntropyLossWeight

% Copyright 2019 The MathWorks Inc.

% policy gradient loss
PolicyGradientLoss = crossentropy(Policy, reshape(LossVariable.Advantage,size(Policy)));

% entropy loss
EntropyLoss = rl.loss.policyEntropyDiscrete(Policy, ...
    LossVariable.EntropyLossWeight,numel(LossVariable.Advantage));

% total loss 
Loss = PolicyGradientLoss + EntropyLoss;
end