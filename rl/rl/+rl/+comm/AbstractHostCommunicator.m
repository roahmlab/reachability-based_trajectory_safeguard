classdef AbstractHostCommunicator < handle
% ABSTRACTHOSTCOMMUNICATOR

% Revised: 10-16-2018
% Copyright 2018 The MathWorks, Inc.
    methods (Abstract)
        % callback syntax:
        %   callback(comm,agent,message)
        registerExperienceReceivedCallback(this,agent,cb)
        registerGradientReceivedCallback(this,agent,cb)
        registerTrainingInfoReceivedCallback(this,agent,cb)
        registerParametersRequestedCallback(this,agent,cb)
        sendContinueTraining(this,val,workerID)
        sendParameters(this,params,workerID)
        sendManuallyStopTraining(this)
    end
end