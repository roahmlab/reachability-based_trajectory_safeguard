function Options = rlRepresentationOptions(varargin)
% rlRepresentationOptions: Creates options for RL agent representations
%
%   opt = rlRepresentationOptions returns the default options for RL
%   representations.
%
%   opt = rlRepresentationOptions('Option1',Value1,'Option2',Value2,...)
%   uses name/value pairs to override the default values for
%   'Option1','Option2',...
%
%   Supported options are:
%
%   LearnRate                     Learning rate for the representation
%   Optimizer                     Optimizer for the representation. The default is 'adam'.
%   OptimizerParameters           Applicable parameters for the optimizer.
%   - Momentum                      Contribution of previous step
%   - Epsilon                       Denominator offset
%   - GradientDecayFactor           Decay rate of gradient moving average
%   - SquaredGradientDecayFactor    Decay rate of squared gradient moving average
%   GradientThreshold             Threshold value for gradient
%   GradientThresholdMethod       The method for gradient threshold
%   L2RegularizationFactor        Factor for L2 regularization
%   UseDevice                     Choose cpu, gpu, or others for representation

% Copyright 2017-2018 The MathWorks Inc.

Options = rl.option.rlRepresentationOptions(varargin{:});

end