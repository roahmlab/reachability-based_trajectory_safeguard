function layer = quadraticLayer(varargin)
% QUADRATICLAYER: Create a quadratic layer for learning
%
% LAYER = quadraticLayer creates a quadratic layer that outputs a
% vector of quadratic monomials given an input vector. For example,
% consider an input vector U = [u1 u2 u3]. The output Y = predict(layer,U)
% will result in Y = [u1*u1 u1*u2 u2*u2 u1*u3 u2*u3 u3*u3].
%
% LAYER = quadraticLayer('Option1',Value1,'Option2',Value2,...) uses
% name/value pairs to override the default values for
% 'Option1','Option2',...
%
% Supported options are:
%
%   Name                 Name of the layer
%   Description          Description of the layer
%
% See also: scalingLayer

% Revised: 12-26-2018
% Copyright 2018 The MathWorks, Inc.

try
    layer = rl.layer.QuadraticLayer(varargin{:});
catch ex
    throwAsCaller(ex);
end

