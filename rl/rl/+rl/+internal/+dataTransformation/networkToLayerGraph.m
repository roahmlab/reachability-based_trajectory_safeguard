function lgraph = networkToLayerGraph(net)
% NETWORKTOLAYERGRAPH   Convert a neural network object (SeriesNetwork, 
% DAGNetwork, dlnetwork) to layerGraph object. The output is used as the
% input to rlLayerRepresentation.
%
%   lgraph = rl.internal.dataTransformation.networkToLayerGraph(net) takes an 
%   input network object net (SeriesNetwork or DAGNetwork) and convert it
%   to layerGraph object lgraph. If the input is an array of layers or a
%   layerGraph, return the input as is.

% Copyright 2019 The MathWorks Inc.

validateattributes(net,{'SeriesNetwork','DAGNetwork','dlnetwork',...
    'nnet.cnn.layer.Layer','nnet.cnn.LayerGraph'}, {'nonempty'})

switch class(net)
    case 'SeriesNetwork'
        lgraph = layerGraph(net.Layers);
    case {'DAGNetwork','dlnetwork'}
        lgraph = layerGraph(net);
    case 'nnet.cnn.layer.Layer'
        % add missing name
        assignedLayerNames = arrayfun(@(x) x.Name,net,'UniformOutput',false);
        ix = cellfun(@(x) isempty(x),assignedLayerNames);
        assignedLayerNames(ix) = matlab.lang.makeUniqueStrings(repmat({'LayerName'},sum(ix),1),assignedLayerNames(~ix));
        for ct = 1:length(ix)
            net(ct).Name = assignedLayerNames{ct};
        end
        lgraph = layerGraph(net);
    case 'nnet.cnn.LayerGraph'
        lgraph = net;
end
end