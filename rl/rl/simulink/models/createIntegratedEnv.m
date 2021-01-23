function [env,agtblk,observationInfo,actionInfo] = createIntegratedEnv(mdl,newmdlname,varargin)
% [ENV,AGENTBLOCK,OBSERVATIONSPECS,ACTIONSPECS] = CREATEINTEGRATEDENV(MODEL,NEWMODELNAME)
%
% Create a new Simulink model NEWMODELNAME that references the Simulink
% model MODEL. A Simulink environment object ENV is returned along
% with the path to the agent block AGENTBLOCK, the observation RL data
% specifications OBSERVATIONSPECS, and the action RL data specification
% ACTIONSPEC.
%
% [ENV,AGENTBLOCK,OBSERVATIONSPECS,ACTIONSPECS] = CREATEINTEGRATEDENV(MODEL,NEWMODELNAME,NAME1,VALUE1,...)
%
% Provide additional name-value pairs to control which ports map to
% the observation, action, reward, and isdone signals, as well as control
% how the RL data specifications are created. The supported name-value pairs
% are:
%
%   ObservationPortName         : Name of the observation outport from the
%                                 referenced MODEL (default = "observation")
%   ActionPortName              : Name of the action inport from the
%                                 referenced MODEL (default = "action")
%   RewardPortName              : Name of the reward outport from the
%                                 referenced MODEL (default = "reward")
%   IsDonePortName              : Name of the isdone outport from the
%                                 referenced MODEL (default = "isdone")
%   ObservationBusElementNames  : The observation element names to generate
%                                 RL data specifications for. The provided
%                                 elements must be leaf elements of the
%                                 observation bus. By default, RL data
%                                 specifications will be created for all
%                                 leaf bus elements. This parameter is not
%                                 applicable for observation outports that
%                                 are not bus signals.
%   ObservationDiscreteElements : Specify a finite set of values for an
%                                 observation specification element using a
%                                 name-value pair cell array. The specified
%                                 set of finite values must be castable to
%                                 the element data type.
%   ActionDiscreteElements      : Specify a finite set of values for an
%                                 action specification element using a
%                                 name-value pair cell array. The specified
%                                 set of finite values must be castable to
%                                 the element data type.
%
% Example 1: Create an environment from a referenced model. Specify the
% ports for each signal.
%
%   env = createIntegratedEnv('my_model','integrated_model',...
%       'ObservationPortName','obs',...
%       'ActionPortName','tau',...
%       'RewardPortName','qrreward',...
%       'IsDonePortName','stopsim');
%
% Example 2: Create an environment from a referenced model. The ports are
% labeled according to the default values. Modify the action specification
% to be bounded between [-2 2].
%
%   [~,agentBlk,observationInfo,actionInfo] = createIntegratedEnv('my_model','integrated_model');
%   actionInfo.UpperLimit =  2;
%   actionInfo.LowerLimit = -2;
%   env = rlSimulinkEnv('integrated_model',agentBlk,observationInfo,actionInfo);
%
%
% Example 3: Create an environment from a referenced model. The ports are
% labeled according to the default values. Specify the actions as discrete
% with elements = [-2 1 0 1 2]
%
%   env = createIntegratedEnv('my_model','integrated_model',...
%       'ActionDiscreteElements',{'action',-2:2});
%
% See also: rlSimulinkEnv, bus2RLSpec, rlNumericSpec, rlFiniteSetSpec

% Revised: 8-29-2018
% Copyright 2017-2018 The MathWorks Inc.

%% parse the inputs

p = inputParser();

strval = @(x) validateattributes(x,{'char','string'},{'scalartext'});
addRequired(p,'Model',strval);
addRequired(p,'NewModelName',strval);

addParameter(p,'ObservationPortName','observation',strval);
addParameter(p,'ActionPortName','action',strval);
addParameter(p,'RewardPortName','reward',strval);
addParameter(p,'IsDonePortName','isdone',strval);
addParameter(p,'ObservationDiscreteElements',{},@(x) localValidateObservationsDiscreteElements(x));
addParameter(p,'ActionDiscreteElements'     ,{},@(x)      localValidateActionsDiscreteElements(x));
addParameter(p,'ObservationBusElementNames' ,{},@(x)              localValidateBusElementNames(x));

parse(p,mdl,newmdlname,varargin{:});

s = p.Results;

obsPort = s.ObservationPortName;
actPort = s.ActionPortName;
rwdPort = s.RewardPortName;
isdPort = s.IsDonePortName;

obsDiscreteElements = s.ObservationDiscreteElements;
actDiscreteElements = s.ActionDiscreteElements     ;
obsBusEls           = s.ObservationBusElementNames ;

% make sure the new model does not exist and the newmodel ~= mdl
if exist(fullfile(pwd,[newmdlname,'.slx']),'file')
    error(message('rl:env:CreateIntegratedEnvNewModelAlreadyExists',newmdlname));
end
if strcmp(newmdlname,mdl)
    error(message('rl:env:CreateIntegratedNewModelSameAsReferencedModel'));
end
% close the model if only loaded
if bdIsLoaded(newmdlname)
    close_system(newmdlname,0);
end

%% make sure the ref model CAN be referenced
numInst = get_param(mdl,'ModelReferenceNumInstancesAllowed');
if strcmp(numInst,'Zero')
    error(message('rl:env:CreateIntegratedReferencedModelNotAllowed',mdl));
end
minAlgLoop = get_param(mdl,'ModelReferenceMinAlgLoopOccurrences');
if strcmp(minAlgLoop,'off')
    warning(message('rl:env:CreateIntegratedMinimizeAlgLoopOff',mdl));
end

%% get the port info
portinfo = rl.util.generateSignalInfoFromModel(mdl,...
    obsPort,actPort,rwdPort,isdPort);

%% create the new model from the template
temppath = fullfile(matlabroot,'toolbox','rl','rl','simulink','models','rlmodeltemplate.slx');
newmdlh = new_system(newmdlname,'FromFile',temppath);
open_system(newmdlh);

envwrapper = [newmdlname,'/Environment'];

envblk = [newmdlname,'/Environment/Environment Model'];
obsblk = [newmdlname,'/Environment/observation'      ];
isdblk = [newmdlname,'/Environment/isdone'           ];
rwdblk = [newmdlname,'/Environment/reward'           ];
actblk = [newmdlname,'/Environment/action'           ];
agtblk = [newmdlname,'/RL Agent'                     ];

obsph = get_param(obsblk,'PortHandles');
actph = get_param(actblk,'PortHandles');
rwdph = get_param(rwdblk,'PortHandles');
isdph = get_param(isdblk,'PortHandles');

% link each model to the same configuration set
cs = getActiveConfigSet(mdl);

cscopy = attachConfigSetCopy(newmdlname,cs,true);
setActiveConfigSet(newmdlname,get_param(cscopy,'Name'));

%% configure the blocks
set_param(envblk,'ModelName',mdl);

%% wire up the signals
envph = get_param(envblk,'porthandles');

% add terminators to unused output ports
usedports = [portinfo.Observation.PortNum,portinfo.Reward.PortNum,portinfo.IsDone.PortNum];
allports = 1:portinfo.NumOutputs;

unusedports = setdiff(allports,usedports);
for i = 1:numel(unusedports)
    port = unusedports(i);
    
    blk = add_block(sprintf('simulink/Sinks/Terminator'),[newmdlname,'/Environment/Terminator'],'MakeNameUnique','on');
    ph = get_param(blk,'porthandles');
    
    snk = ph.Inport(1);
    src = envph.Outport(port);
    
    portpos = get_param(src,'position');
    set_param(blk,'Position',[portpos(1)+20 portpos(2)-15 portpos(1)+50 portpos(2)+15]);
    
    localAddLine(envwrapper,src,snk,'');
end

% observation
src = envph.Outport(portinfo.Observation.PortNum);
snk = obsph.Inport(1);
localAddLine(envwrapper,src,snk,portinfo.Observation.Name);

% reward
src = envph.Outport(portinfo.Reward.PortNum);
snk = rwdph.Inport(1);
localAddLine(envwrapper,src,snk,portinfo.Reward.Name);

% isdone
src = envph.Outport(portinfo.IsDone.PortNum);
snk = isdph.Inport(1);
localAddLine(envwrapper,src,snk,portinfo.IsDone.Name);

% action
snk = envph.Inport(portinfo.Action.PortNum);
src = actph.Outport(1);
localAddLine(envwrapper,src,snk,portinfo.Action.Name);

%% create the observtion/action specs and the environment
try
    % try to beutify the model
    % glm = Simulink.internal.variantlayout.LayoutManager(newmdlname);
    glm = Simulink.internal.variantlayout.LayoutManager(envwrapper);
    layoutModel(glm);
    
    [observationInfo,actionInfo] = localGenerateSpecDataFromPortInfo(portinfo,newmdlname,obsDiscreteElements,actDiscreteElements,obsBusEls);
    env = rlSimulinkEnv(newmdlname,agtblk,observationInfo,actionInfo);
    % flush the event queue
    drawnow();
catch ex
    close_system(newmdlname,0);
    rethrow(ex);
end


function [oinfo,ainfo] = localGenerateSpecDataFromPortInfo(portinfo,mdl,odiscreteElements,adiscreteElements,obsBusEls)
oinfo = localGenerateSpecDataFromStruct(portinfo.Observation,mdl,odiscreteElements,obsBusEls);
ainfo = localGenerateSpecDataFromStruct(portinfo.Action     ,mdl,adiscreteElements,{});

function info = localGenerateSpecDataFromStruct(s,mdl,discreteElements,busEls)
try
    if s.IsBus
        info = bus2RLSpec(s.DataType,...
            'Model',mdl,...
            'DiscreteElements',discreteElements,...
            'BusElementNames',busEls);
    else
        if isempty(discreteElements)
            info = rlNumericSpec(s.Dimensions);
        else
            if ~strcmp(discreteElements{1},s.Name)
                error(message('rl:env:CreateIntegratedEnvNoMatchingPortSpec',discreteElements{1},s.Name));
            end
            info = rlFiniteSetSpec(cast(discreteElements{2},s.DataType));
        end
        info.Name = s.Name;
    end
catch ex
    me = MException(message('rl:env:CreateIntegratedEnvErrorCreatingSpec',s.Name));
    me = addCause(me,ex);
    throwAsCaller(me);
end

function localAddLine(mdl,src,snk,name)
line = add_line(mdl,src,snk,'autorouting','smart');
set_param(line,'Name',name);

function localValidateObservationsDiscreteElements(val)
if ~isempty(val)
    if ~iscell(val)
        error(message('rl:env:CreateIntegratedEnvInvalidObservationDiscreteElementsArg'))
    end
    N = numel(val);
    if mod(N,2)
        error(message('rl:env:CreateIntegratedEnvInvalidObservationDiscreteElementsArg'));
    end
    names = val(1:2:end-1);
    if ~(isstring(names) || iscellstr(names))
        error(message('rl:env:CreateIntegratedEnvInvalidObservationDiscreteElementsArg'));
    end
end

function localValidateActionsDiscreteElements(val)
if ~isempty(val)
    if ~iscell(val)
        error(message('rl:env:CreateIntegratedEnvInvalidActionDiscreteElementsArg'))
    end
    N = numel(val);
    if mod(N,2)
        error(message('rl:env:CreateIntegratedEnvInvalidActionDiscreteElementsArg'));
    end
    names = val(1:2:end-1);
    if ~(isstring(names) || iscellstr(names))
        error(message('rl:env:CreateIntegratedEnvInvalidActionDiscreteElementsArg'));
    end
end

function localValidateBusElementNames(val)
if ~isempty(val)
    if isstring(val)
        val = cellstr(val);
    end
    if ~iscellstr(val)
        error(message('rl:env:CreateIntegratedEnvInvalidObservationBusElementNamesArg'))
    end
end

