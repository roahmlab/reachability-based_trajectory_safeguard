function xBounded = boundAwayFromZero(x)
% boundAwayFromZero   Round any values less eps to eps
%
%   xBounded = boundAwayFromZero(x) takes an input array x and returns an
%   array xBounded where all values less than eps are replaced with eps.

% Copyright 2019 The MathWorks Inc.

xBounded = x;
precision = 'single'; % assume always single
xBounded(abs(xBounded) < eps(precision)) = eps(precision);
end