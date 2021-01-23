classdef Task < handle
% TASK
%
% A Task is a base class for running code on a thread via the run method.
% TaskSpecs will build Tasks on threads/workers to be run

% Revised: 6-13-2019
% Copyright 2019 The MathWorks, Inc.

    properties (Access = protected)
        WorkerID
    end
    methods 
        function this = Task()
            this.WorkerID = rl.util.getWorkerID();
        end
    end
    methods (Sealed)
        function varargout = run(this)
            % run the task
            [varargout{1:nargout}] = runImpl(this);
        end
    end
    methods (Access = protected,Abstract)
        varargout = runImpl(this)
    end
end