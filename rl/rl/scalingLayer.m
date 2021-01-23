function layer = scalingLayer(varargin)
% SCALINGLAYER: Create a scaling layer for learning
%
% LAYER = scalingLayer creates a scaling layer that linearly scales and
% biases an input array U so the output Y = scale.*U + bias. This layer is
% useful for scaling and shifting outputs of nonlinear layers, such as tanh
% and sigmoid. By default, the scale is 1 and the bias is 0.
%
% LAYER = scalingLayer('Option1',Value1,'Option2',Value2,...) uses
% name/value pairs to override the default values for
% 'Option1','Option2',...
%
% Supported options are:
%
%   Name                 Name of the layer
%   Description          Description of the layer
%   Scale                Elementwise scale of the input (default = 1)
%   Bias                 Elementwise bias of the input (default = 0)
%
% See also: quadraticLayer

% Revised: 12-26-2018
% Copyright 2018 The MathWorks, Inc.

try
    layer = rl.layer.ScalingLayer(varargin{:});
catch ex
    throwAsCaller(ex);
end
