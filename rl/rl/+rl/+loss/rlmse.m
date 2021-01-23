function Loss = rlmse(ModelOutput,Target)
% Mean squared error loss.
%   ModelOutput: dlarray, output from model
%   Target: mse target, can be numeric

% Copyright 2019 The MathWorks Inc.

Loss = mse(ModelOutput, reshape(Target,size(ModelOutput)));

end