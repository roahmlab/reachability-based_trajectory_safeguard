function str = handleMultipleOutputStrings(str)
% HANDLEMULTIPLEOUTPUTSTRINGS
% u1,u2,u3 -> [u1,u2,u3]

% Revised: 2-4-2019
% Copyright 2019 The MathWorks, Inc.
if ismember(',',str)
    str = ['[',str,']'];
end
