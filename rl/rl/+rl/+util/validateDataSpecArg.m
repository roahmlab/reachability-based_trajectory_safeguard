function validateDataSpecArg(DataSpec, Filename, VarName, ArgIndex)
%VALIDATEDATASPECARG Type, size checking for representation constructor
% data spec

% Copyright 2019 The MathWorks, Inc.

validateattributes(DataSpec, {'rl.util.RLDataSpec'}, {'vector', 'nonempty'}, Filename, VarName, ArgIndex);
end