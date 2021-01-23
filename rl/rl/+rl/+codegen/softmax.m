function p = softmax(n)
% SOFTMAX
% map the given vector n to probabilities vector p sums to 1,
% with elements exponentially proportional to the respective elements in n.

% Copyright 2019 The MathWorks, Inc.

%#codegen

v = exp(n - max(n));
p = v/sum(v);
