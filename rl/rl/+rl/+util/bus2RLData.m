function data = bus2RLData(busName,varargin)
% BUS2RLDATA

% Revised: 10-27-2018
% Copyright 2018 The MathWorks, Inc.

% parse args
p = inputParser;
strval = @(x) validateattributes(x,{'char','string'},{'scalartext'});
addRequired(p,'busName',strval);
addParameter(p,'Model','',strval);
addParameter(p,'BusElementNames',{},@(x) localValidateBusElementNames(x));
addParameter(p,'DiscreteElements',{},@(x) localValidateDiscreteElements(x));

parse(p,busName,varargin{:});

r = p.Results;

mdl              = r.Model           ;
busElementNames  = r.BusElementNames ;
discreteElements = r.DiscreteElements;

denames = discreteElements(1:2:end-1);
devals  = discreteElements(2:2:end  );

% get the bus object
try
    if isempty(mdl)
        busObj = evalin('base',busName);
    else
        busObj = Simulink.data.evalinGlobal(mdl,busName);
    end
catch ex
    if isempty(mdl)
        me = MException(message('rl:general:Bus2RLDataCouldNotFindBusInBase',busName));
    else
        me = MException(message('rl:general:Bus2RLDataCouldNotFindBusInGlobal',busName));
    end
    me = addCause(me,ex);
    throw(me);
end
elements = getLeafBusElements(busObj);

% get the element names, data types, dimensions, min/max, and description
elnames = {elements.Name       }';
dts     = {elements.DataType   }';
dims    = {elements.Dimensions }';
mins    = {elements.Min        }';
maxs    = {elements.Max        }';
descs   = {elements.Description}';

% make sure the provided element names match the element names
if isempty(busElementNames)
    busElementNames = elnames;
end
elidx = ismember(busElementNames,elnames);
if ~all(elidx)
    f1 = find(~elidx,1);
    error(message('rl:general:Bus2RLDataElementDoesNotExist',busElementNames{f1},busName));
end

elidx = find(ismember(elnames,busElementNames));
N = numel(elidx);
discreteElements = cell(N,1);
if ~isempty(denames)
    for i = 1:numel(denames)
        den = denames{i};
        deidx = ismember(elnames,den);
        if ~any(deidx)
            error(message('rl:general:Bus2RLDataElementDoesNotExist',den,busName));
        end
        discreteElements{deidx} = devals{i};
    end
end
if numel(discreteElements) ~= N
    error(message('rl:general:Bus2RLDataDiscreteElementsMismatch',numel(discreteElements),N));
end

% create the rl data for each element
for i = 1:N
    idx = elidx(i);
    de = discreteElements{i};
    
    % extract
    elname  = elnames{idx};
    dt      = dts    {idx};
    dim     = dims   {idx};
    desc    = descs  {idx};
    min_    = localConvertEmptyToInf(mins{idx},dim);
    max_    = localConvertEmptyToInf(maxs{idx},dim);
    
    ii = isinf(min_);
    min_(ii) = -min_(ii);
    
    try
        if isempty(de)
            % create NDim data
            d = rlNumericSpec(dim,...
                'DataType'  ,dt,...
                'LowerLimit',min_,....
                'UpperLimit',max_);
        else
            % create discrete data
            d = rlFiniteSetSpec(cast(de,dt));
        end
    catch ex
        me = MException(message('rl:general:Bus2RLDataErrorCreatingSpecData',elname));
        me = addCause(me,ex);
        throw(me);
    end
    
    d.BusName     = busName ;
    d.Name        = elname  ;
    d.Description = desc    ;
    
    data(i) = d; %#ok<AGROW>
end

function val = localConvertEmptyToInf(val,dim)
if isempty(val)
    val = Inf.*ones(dim);
end

function localValidateBusElementNames(val)
if ~isempty(val)
    if isstring(val)
        val = cellstr(val);
    end
    if ~iscellstr(val)
        error(message('rl:general:Bus2RLDataInvalidBusElementNamesArg'))
    end
end

function localValidateDiscreteElements(val)
if ~isempty(val)
    if ~iscell(val)
        error(message('rl:general:Bus2RLDataInvalidDiscreteElementsArg'))
    end
    N = numel(val);
    if mod(N,2)
        error(message('rl:general:Bus2RLDataInvalidDiscreteElementsArg'));
    end
    names = val(1:2:end-1);
    if ~(isstring(names) || iscellstr(names))
        error(message('rl:general:Bus2RLDataInvalidDiscreteElementsArg'));
    end
end


