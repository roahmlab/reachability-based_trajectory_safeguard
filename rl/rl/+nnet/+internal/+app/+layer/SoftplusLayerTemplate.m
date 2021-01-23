classdef SoftplusLayerTemplate < nnet.internal.app.plugins.LayerTemplate
% SOFTPLUSLAYERTEMPLATE Represents a Softplus Layer in DeepNetworkDesigner app
%

% Copyright 2019 The MathWorks, Inc

    properties
        ConstructorName = "softplusLayer"
        RequiredArguments = []
        OptionalArguments = struct('Name', 'softplus')
    end
    methods
        function configure(this)
            this.Group = nnet.internal.app.plugins.layer.LayerGroup.Activation;
            this.IconPath = iIconPath();
        end
    end
    
end

function iconPath = iIconPath()
iconPath = fullfile("toolbox", "nnet", "deepapp", "web", "editor", ...
    "deepapp-editor", "images", "palette", "rl-layer-SoftplusLayer.svg");
end