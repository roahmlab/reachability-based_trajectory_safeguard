function portinfo = generateSignalInfoFromModel(mdl,...
    obsPort,actPort,rwdPort,isdPort)
% GENERATESIGNALINFOFROMMODEL

% Revised: 12-4-2018
% Copyright 2018 The MathWorks, Inc.

%% compile the model and get port data
wascompiled = strcmp(get_param(mdl,'SimulationStatus'),'paused');
clnup = onCleanup(@()localCleanup(mdl,wascompiled));
if ~wascompiled
    feval(mdl,[],[],[],'compile');
end

rootU = find_system(mdl,'SearchDepth',1,'BlockType','Inport' );
rootY = find_system(mdl,'SearchDepth',1,'BlockType','Outport');

nu = numel(rootU);
ny = numel(rootY);

nameU = get_param(rootU,'Name');
nameY = get_param(rootY,'Name');

dtU = cellfun(@(h) get_param(h.Outport,'CompiledPortDataType'),get_param(rootU,'porthandles'),'UniformOutput',false);
dtY = cellfun(@(h) get_param(h.Inport ,'CompiledPortDataType'),get_param(rootY,'porthandles'),'UniformOutput',false);

dimU = cellfun(@(h) get_param(h.Outport,'CompiledPortDimensions'),get_param(rootU,'porthandles'),'UniformOutput',false);
dimY = cellfun(@(h) get_param(h.Inport ,'CompiledPortDimensions'),get_param(rootY,'porthandles'),'UniformOutput',false);

isbusU = ismember(cellfun(@(h) get_param(h.Outport,'CompiledBusType'),get_param(rootU,'porthandles'),'UniformOutput',false),'NON_VIRTUAL_BUS');
isbusY = ismember(cellfun(@(h) get_param(h.Inport ,'CompiledBusType'),get_param(rootY,'porthandles'),'UniformOutput',false),'NON_VIRTUAL_BUS');

%% build port info

% get the idx to the appropriate port
obsIdx = find(ismember(nameY,obsPort)); 
actIdx = find(ismember(nameU,actPort));
rwdIdx = find(ismember(nameY,rwdPort));
isdIdx = find(ismember(nameY,isdPort));

% error out if we can't find matching ports
if isempty(obsIdx)
    error(message('rl:env:NoMatchingOutport',obsPort,mdl));
end
if isempty(actIdx)
    error(message('rl:env:NoMatchingInport',actPort,mdl));
end
if isempty(rwdIdx)
    error(message('rl:env:NoMatchingOutport',rwdPort,mdl));
end
if isempty(isdIdx)
    error(message('rl:env:NoMatchingOutport',isdPort,mdl));
end

portinfo.NumInputs = nu;
portinfo.NumOutputs = ny;

portinfo.Observation.PortNum = obsIdx;
portinfo.Observation.DataType = dtY{obsIdx};
portinfo.Observation.Dimensions = localGetSize(dimY{obsIdx});
portinfo.Observation.Name = obsPort;
portinfo.Observation.IsBus = isbusY(obsIdx);

portinfo.Action.PortNum = actIdx;
portinfo.Action.DataType = dtU{actIdx};
portinfo.Action.Dimensions = localGetSize(dimU{actIdx});
portinfo.Action.Name = actPort;
portinfo.Action.IsBus = isbusU(actIdx);

portinfo.Reward.PortNum = rwdIdx;
portinfo.Reward.DataType = dtY{rwdIdx};
portinfo.Reward.Dimensions = localGetSize(dimY{rwdIdx});
portinfo.Reward.Name = rwdPort;
portinfo.Reward.IsBus = false;

portinfo.IsDone.PortNum = isdIdx;
portinfo.IsDone.DataType = dtY{isdIdx};
portinfo.IsDone.Dimensions = localGetSize(dimY{isdIdx});
portinfo.IsDone.Name = isdPort;
portinfo.IsDone.IsBus = false;

localCheckForScalarSignal('reward',portinfo.Reward.Dimensions);
localCheckForScalarSignal('isdone',portinfo.IsDone.Dimensions);

function localCleanup(mdl,wascompiled)
iscompiled = strcmp(get_param(mdl,'SimulationStatus'),'paused');
if iscompiled && ~wascompiled
    feval(mdl,[],[],[],'term');
end

function sz = localGetSize(dim)
if localIsBus(dim)
    sz = [1 1];
else
    sz = dim(2:end);
    if numel(sz) < 2
        sz = [sz 1];
    end
end

function localCheckForScalarSignal(type,dim)
if localIsBus(dim) || prod(dim(2:end)) > 1
    error(message('rl:env:SignalMustBeScalarNonBus',type));
end

function val = localIsBus(dim)
% see https://www.mathworks.com/help/simulink/ug/determining-output-signal-dimensions.html#bsud_ek-10
val = dim(1) == -2;

