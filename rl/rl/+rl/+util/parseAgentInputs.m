function [rep,opt,UseDefault] = parseAgentInputs(Arguments,OptionsClass,DefaultOptions)
% validate RL agent inputs are representations and valid options

% Copyright 2018 The MathWorks, Inc.

% representation check
rep = Arguments(cellfun(@(x) isa(x,'rl.util.rlAbstractRepresentation'),Arguments));
if isempty(rep)
    error(message('rl:agent:errRequiredRepresentationInput'))    
end

% options check
UseDefault = false;
opt = Arguments(cellfun(@(x) isa(x,'rl.option.AgentGeneric'),Arguments));

if numel(Arguments)~=( numel(rep)+numel(opt) )
    error(message('rl:agent:errInvalidAgentInput'));
end

if isempty(opt)
    opt{1} = DefaultOptions;
    UseDefault = true;
else
    % check option is compatible
    if ~isa(opt{1},OptionsClass)
        error(message('rl:agent:errMismatchedOption'));        
    end
end