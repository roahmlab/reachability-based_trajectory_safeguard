function EntropyLoss = policyEntropyDiscrete(Probability, EntropyLossWeight, NumExperience)
% Discrete policy entropy loss to encourage policy exploration
% Entropy = -sum(Prob(x).*log(Prob(x)))
% EntropyLoss = EntropyLossWeight * sum(Entropy)/NumObs

% Copyright 2019 The MathWorks Inc.

if EntropyLossWeight
    EntropyLoss = EntropyLossWeight * ...
        sum(Probability .* log(rl.internal.dataTransformation.boundAwayFromZero(Probability)),'all')./NumExperience;
else
    EntropyLoss = 0;
end

end