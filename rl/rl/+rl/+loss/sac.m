function Loss = sac(ProbabilityParams, LossVariable)
% Soft Actor Critic Policy Loss.
% ProbabilityParams: mean and variance from the stochastic actor representation
% LossVariable: Includes all the necessary objects/variables for
% calculating the loss.
%
% Copyright 2019 The MathWorks, Inc.

Actions = sample(LossVariable.Actor.SamplingStrategy, {ProbabilityParams});
SquashedActions = tanh(Actions{1});
AdjustedActions = SquashedActions * LossVariable.ActionScale + LossVariable.ActionShift;
AdjustedActions = {AdjustedActions};
Densities = evaluate(LossVariable.Actor.SamplingStrategy, ProbabilityParams, Actions{1});
QValue = getValue(LossVariable.Critic1, LossVariable.Observations, AdjustedActions);
if ~isempty(LossVariable.Critic2)
    QValue2 = getValue(LossVariable.Critic2, LossVariable.Observations, AdjustedActions);
    QValue = min(QValue, QValue2);
end
Densities = prod(Densities, 1);
LogDensities = log(rl.internal.dataTransformation.boundAwayFromZero(Densities));
PDFCorrection = sum(log(rl.internal.dataTransformation.boundAwayFromZero((1 - tanh(Actions{1}).^2) .* LossVariable.ActionScale)), 1);
LogDensities = LogDensities - PDFCorrection;
Loss = LossVariable.EntropyWeight * LogDensities - QValue;
Loss = mean(Loss, 'all');
end