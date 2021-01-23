function idx = getElementIndex(set,data)
% GETELEMENTINDEX
% given a set and data, determine the index of the data. Will
% return [] if data cannot be found in the set

% Revised: 2-8-2019
% Copyright 2019 The MathWorks, Inc.

%#codegen

if iscell(set)
    idx = [];
    for i = 1:numel(set)
        if isequal(set{i},data)
            idx = i;
            break;
        end
    end
else
    idx = find(ismember(set,data),1);
end

