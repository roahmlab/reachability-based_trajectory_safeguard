function Model = createInternalModelFactory(Model, Options, ObservationNames, ActionNames, InputSize, OutputSize)
% CREATEINTERNALMODELFACTORY Create an internal model from user-defined model

% Copyright 2019 The MathWorks, Inc.

% InputSize and OutputSize only applicable to basis function model

UseDevice = Options.UseDevice;

if ~isa(Model, 'rl.representation.model.rlAbstractModel')
    switch class(Model)
        case {'nnet.cnn.layer.Layer','nnet.cnn.LayerGraph','SeriesNetwork','DAGNetwork','dlnetwork'}
            Model = rl.internal.dataTransformation.networkToLayerGraph(Model);
            Model = rl.representation.model.rlLayerModel(Model, UseDevice, ObservationNames, ActionNames);
        case {'rlTable'}
            Model = rl.representation.model.rlTableModel(Model, UseDevice);
        case {'cell'}
            Model = rl.representation.model.rlBasisFunctionModel(Model, UseDevice, InputSize, OutputSize);
        otherwise
            error(message('rl:agent:errUnknownModel'));
    end
end
end