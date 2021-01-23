classdef AbstractWorkerCommunicator < handle
% ABSTRACTWORKERCOMMUNICATOR

% Revised: 10-16-2018
% Copyright 2018 The MathWorks, Inc.

    methods (Abstract)
        params = handshakeExperience(this,exp)
        params = handshakeGradient(this,gradientBuffer)
        params = handshakeParameterRequest(this)
        params = receiveParameters(this)
        continueTraining = handshakeTrainingInfo(this,info)
        sendTrainingInfo(this,info)
        cont = receiveContinueTraining(this,varargin)
        sendExperience(this,exp)
        sendGradient(this,gradientBuffer)
        id = getID(this)
        registerParametersReceivedCallback(this,agent,cb)
        registerContinueTrainingCallback(this,agent,cb)
        registerManuallyStoppedTrainingCallback(this,agent,cb)
    end
end