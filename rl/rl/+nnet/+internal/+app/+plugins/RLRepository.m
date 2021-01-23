classdef RLRepository < nnet.internal.app.plugins.Repository
% RLREPOSITORY Defines content from RL toolbox to include in DeepNetworkDesigner app
%

% Copyright 2019 The MathWorks, Inc

    properties
        Layers = {...
            % Activation
            nnet.internal.app.layer.SoftplusLayerTemplate();
            
            % Normalization and Utility
            nnet.internal.app.layer.ScalingLayerTemplate();
            nnet.internal.app.layer.QuadraticLayerTemplate();
            }
    end
end