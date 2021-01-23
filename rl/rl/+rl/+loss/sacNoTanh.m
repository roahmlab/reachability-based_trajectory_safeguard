function Loss = sacNoTanh(ProbabilityParams, LossVariable)
% Soft Actor Critic Policy Loss.
% ProbabilityParams: mean and variance from the stochastic actor representation
% LossVariable: Includes all the necessary objects/variables for
% calculating the loss.
%
% Copyright 2019 The MathWorks, Inc.

Actions = sample(LossVariable.Actor.SamplingStrategy, {ProbabilityParams});
AdjustedActions = LossVariable.ActionInfo.saturate(Actions);
Densities = evaluate(LossVariable.Actor.SamplingStrategy, ProbabilityParams, AdjustedActions{1});
QValue = getValue(LossVariable.Critic1, LossVariable.Observations, AdjustedActions);
if ~isempty(LossVariable.Critic2)
    QValue2 = getValue(LossVariable.Critic2, LossVariable.Observations, AdjustedActions);
    QValue = min(QValue, QValue2);
end
Densities = prod(Densities, 1);
LogDensities = log(rl.internal.dataTransformation.boundAwayFromZero(Densities));
Loss = LossVariable.EntropyWeight * LogDensities - QValue;
Loss = mean(Loss, 'all');
end