function [Option,ObservationNames,ActionNames] = parseRepresentationInput(varargin)
% Parse representation option, observation and action name if any

% Copyright 2019 The MathWorks, Inc.

% NOTE: manual parsing since inputs are not all name-value pairs
ix = find(cellfun(@(x) isa(x,'rl.option.rlRepresentationOptions'),varargin),1);
if isempty(ix)
    Option = rlRepresentationOptions();
else
    Option = varargin{ix};
end

% parse observation name
ObservationNames = {};
ix = find(cellfun(@(x) (ischar(x) || isstring(x)) && all(strcmpi(x,'Observation')),varargin),1);
if ~isempty(ix) && ix < length(varargin)
    ObservationNames = varargin{ix+1};
end
if ~isempty(ObservationNames)
    ObservationNames = iStandardizeName(ObservationNames);
end

% parse action name
ActionNames = {};
ix = find(cellfun(@(x) (ischar(x) || isstring(x)) && all(strcmpi(x,'Action')),varargin),1);
if ~isempty(ix) && ix < length(varargin)
    ActionNames = varargin{ix+1};
end
if ~isempty(ActionNames)
    ActionNames = iStandardizeName(ActionNames);
end
end

function Names = iStandardizeName(Names)
iCheckValidNameInput(Names);
switch class(Names)
    case 'string'
        Names = cellstr(Names);
    case 'char'
        Names = {Names};
    case 'cell'
        Names = cellfun(@char,Names,'UniformOutput',false);
end
Names = reshape(Names,1,[]);
end

function iCheckValidNameInput(Names)
% cell array of char array; char array; string array
if iscell(Names)
    IsValidNameInput = all(cellfun(@(x) ischar(x), Names));
else
    IsValidNameInput = ischar(Names) || isstring(Names);
end
if ~IsValidNameInput
    error(message('rl:agent:errInvalidNeuralNetIOName'));
end
end