function spec = rlNumericSpec(Dimension, varargin)
% rlNumericSpec: Creates an object to define the action or observation space.
%
%   spec = rlNumericSpec(Dimension) creates data spec with the shape
%   defined by the vector Diminesion.
%
%   For example:
%     spec = rlNumericSpec([160,210,3],'LowerLimit', 0, 'UpperLimit',1);
%
%   Use "usample" to create a sample from this specification. For example:
%     s = usample(spec)
%
%   See also: rlFiniteSetSpec

% Copyright 2018 The MathWorks Inc.

spec = rl.util.rlNumericSpec(Dimension, varargin{:});

end