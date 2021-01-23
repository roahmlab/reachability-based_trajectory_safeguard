function id = getWorkerID(hasPCT)
% GETWORKERID

% Revised: 3-25-2019
% Copyright 2019 The MathWorks, Inc.

% get an id for workers. If parallel is not running id = 0

if ~nargin
    hasPCT = matlab.internal.parallel.isPCTInstalled();
end

id = 0;
if hasPCT
    task = getCurrentTask();
    if ~isempty(task)
        id = task.ID;
    end
end

