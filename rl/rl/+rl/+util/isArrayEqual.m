function flag = isArrayEqual(a1,a2)
% Helper function to compare two arrays with different sizes of first three
% dimensions.
% 
% flag = rl.util.isArrayEqual([2 1],[2 1 1])

% Copyright 2018-2019 The MathWorks Inc.

% get sizes
n1=length(a1);
n2=length(a2);

% make them equal
a1 = [a1 ones(1,3-n1)];
a2 = [a2 ones(1,3-n2)];

flag = isequal(a1(1:3),a2(1:3));

end
