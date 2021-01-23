classdef RLEventData < event.EventData
% EXPERIENCERECEVIEDEVENTDATA

% Copyright 2018 The MathWorks, Inc.

    properties (SetAccess = immutable)
        Data
    end
    methods
        function this = RLEventData(data)
            this.Data = data;
        end
    end
end