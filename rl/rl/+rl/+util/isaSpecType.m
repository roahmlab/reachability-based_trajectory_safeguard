function Flag = isaSpecType(DataSpec, SpaceIdentifier)
%ISASPECTYPE returns true if all DataSpec (can be scalar, vector or matrix)
%is the same type as SpaceIdentifier. SpaceIdentifier can be 'discrete' or
%'continuous'.
% Example: 
%   spec1 = rlFiniteSetSpec([1 1]);
%   spec2 = rlFiniteSetSpec([1;2;3;4]);
%   spec = [spec1 spec2];
%   isaSpecType(spec,'discrete') % return true

% Copyright 2019 The MathWorks, Inc.

switch SpaceIdentifier
    case 'discrete'
        ClassName = 'rl.util.rlFiniteSetSpec';
    case 'continuous'
        ClassName = 'rl.util.rlNumericSpec';
    otherwise
        error(message('rl:general:errInvalidSpaceType'));
end
DataSpec = DataSpec(:);
Flag = all(arrayfun(@(x) isa(x,ClassName),DataSpec));

end