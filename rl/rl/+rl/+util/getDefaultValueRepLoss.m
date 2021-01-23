function Loss = getDefaultValueRepLoss()
%GETDEFAULTVALUEREPLOSS Returns the default loss function handle for value 
%representation
%   Default loss is mean squared error

% Copyright 2019 The MathWorks, Inc.

Loss = @rl.loss.rlmse;
end

