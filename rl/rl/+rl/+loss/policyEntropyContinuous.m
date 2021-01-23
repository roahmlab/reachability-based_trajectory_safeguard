function EntropyLoss = policyEntropyContinuous(MeanAndStd, EntropyLossWeight, NumExperience)
% Differential policy entropy loss to encourage policy exploration
% Entropy = -0.5*log(2*pi*exp*Std^2)
% EntropyLoss = EntropyLossWeight * sum(Entropy)/NumObs

% Copyright 2019 The MathWorks Inc.

% collect variance
ParamDim = finddim(MeanAndStd,'C');
NumAction = size(MeanAndStd, ParamDim)/2;
Std = rl.internal.dataTransformation.generalSubref(MeanAndStd,(NumAction+1:NumAction*2),ParamDim);

% compute entropy loss
Entropy = -0.5*log(rl.internal.dataTransformation.boundAwayFromZero(2*pi*exp(1)*Std.^2));
EntropyLoss = EntropyLossWeight * sum(Entropy,'all')./NumExperience;

end