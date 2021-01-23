function layer = softplusLayer(varargin)
% SOFTPLUSLAYER     Softplus layer
%
%   LAYER = softplusLayer() creates a softplus layer that implements 
%   softplus activation Y = log(1+exp(X)) to ensure the output is always 
%   positive.
%   This layer is useful when constructing continuous Gaussian policy deep 
%   neural networks, for which the standard deviation output must be positive.
%
%   LAYER = softplusLayer('Option1',Value1,'Option2',Value2) specifies 
%   optional name/value pairs for creating the layer:
%
%       'Name'                 - Name of the layer
%       'Description'          - Description of the layer
%
%   Example of continuous gaussian policy deep neural network:
%   net = [
%       imageInputLayer([4 1 1], 'Normalization', 'none', 'Name', 'state')
%       fullyConnectedLayer(24, 'Name', 'CommonFC')
%       reluLayer('Name', 'CommonRelu')
%       ];
%   meanPath = [
%       fullyConnectedLayer(1, 'Name', 'meanFC')
%       tanhLayer('Name','tanh')
%       scalingLayer('Name','mean','Scale',10)
%       ];
%   stdPath = [
%       fullyConnectedLayer(1, 'Name', 'stdFC')
%       softplusLayer('Name', 'std')
%       ];
%   concatPath = concatenationLayer(3,2,'Name','actionProbParam');
% 
%   net = layerGraph(net);
%   net = addLayers(net, meanPath);
%   net = addLayers(net, stdPath);
%   net = addLayers(net, concatPath);
%   net = connectLayers(net,'CommonRelu','meanFC/in');
%   net = connectLayers(net,'CommonRelu','stdFC/in');
%   net = connectLayers(net,'mean','actionProbParam/in1');
%   net = connectLayers(net,'std','actionProbParam/in2');
%
% See also: quadraticLayer, scalingLayer.

% Copyright 2019 The MathWorks, Inc.

try
    layer = rl.layer.SoftplusLayer(varargin{:});
catch ex
    throwAsCaller(ex);
end
