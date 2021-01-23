classdef ScalingLayerTemplate < nnet.internal.app.plugins.LayerTemplate
% SCALINGLAYERTEMPLATE Represents a Scaling Layer in DeepNetworkDesigner app
%

% Copyright 2019 The MathWorks, Inc

    properties
        ConstructorName = "scalingLayer"
        RequiredArguments = []
        OptionalArguments = struct('Name', 'scaling')
    end
    methods
        function configure(this)
            this.Group = nnet.internal.app.plugins.layer.LayerGroup.NormalizationAndUtility;
            this.IconPath = iIconPath();
        end
    end
end

function iconPath = iIconPath()
iconPath = fullfile("toolbox", "nnet", "deepapp", "web", "editor", ...
    "deepapp-editor", "images", "palette", "rl-layer-ScalingLayer.svg");
end