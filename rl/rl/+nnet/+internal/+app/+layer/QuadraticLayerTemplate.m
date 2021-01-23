classdef QuadraticLayerTemplate < nnet.internal.app.plugins.LayerTemplate
% QUADRATICLAYERTEMPLATE Represents a Quadratic Layer in DeepNetworkDesigner app
%

% Copyright 2019 The MathWorks, Inc

    properties
        ConstructorName = "quadraticLayer"
        RequiredArguments = []
        OptionalArguments = struct('Name', 'quadratic')
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
    "deepapp-editor", "images", "palette", "rl-layer-QuadraticLayer.svg");
end