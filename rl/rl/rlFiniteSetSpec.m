function spec = rlFiniteSetSpec(Elements)
% rlFiniteSetSpec: Creates an object to define the finite set of actions or observations.
%
%   spec = rlFiniteSetSpec(Elements) creates data spec with the discrete set 
%   defined by the elements Elements.
%   
%   For example:
%     spec = rlFiniteSetSpec([0;1;2;3])
%     spec = rlFiniteSetSpec({[0,1];[1,1];[1,2];[1,3]})
%
%   Use "usample" to create a sample from this specification. For example:
%      s = usample(spec)
%
%   See also: rlNumericSpec

% Copyright 2018 The MathWorks Inc.

spec = rl.util.rlFiniteSetSpec(Elements);

end