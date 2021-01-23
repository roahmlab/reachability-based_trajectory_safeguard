function opt = EntropyWeightOptions(varargin)
% EntropyWeightOptions: Creates options for Entropy Weight
%
%   OPT = EntropyWeightOptions returns the default options for
%   EntropyWeight
%    
%   OPT = EntropyWeightOptions('Option1',Value1,'Option2',Value2,...) uses
%   name/value pairs to override the default values for 'Options1',
%   'Option2', ...
%
%   Supported options are:
%
%   EntropyWeight                       Weight of the entropy component
%   TargetEntropy                       Desired Target Entropy of the policy
%   Optimizer                           Name of the entropy optimizer
%   LearnRate                           The learning rate for the optimizer
%   GradientThreshold                   Value of threshold for the gradient
%   OptimizerParameters                 Optimizer specific parameters
opt = rl.option.EntropyWeightOptions(varargin{:});
end