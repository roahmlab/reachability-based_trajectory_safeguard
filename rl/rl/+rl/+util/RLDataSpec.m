classdef (Abstract = true) RLDataSpec < matlab.mixin.Heterogeneous
    % RLData Abstract data class for action and observation space
    % definitions 
    
    % Copyright 2018-2019 The MathWorks Inc.
    
    %%
    properties
        Name string
        Description string
    end
    
    %%
    properties (SetAccess = protected)
        % Dimension of the data space
        Dimension 
        
        % Data type of the data space 
        DataType (1,1 )string  {mustBeMember(DataType, ["double", "single", "int8","int16","int32","int64","uint8","uint16","uint32","uint64"])} =  "double"
    end
    
    %%
    properties (Hidden)
        BusName string
    end
    
    %%
    methods (Access = protected)
        function obj = RLDataSpec(Dimension, DataType)
            % RLData constructor
            obj.Dimension = Dimension;
            if nargin == 2
                obj.DataType = DataType;
            end
        end
        
    end
    
    %%
    methods (Hidden,Sealed)
        function [val,busName] = isBusSpec(this)
            % check to see if the spec is sourced from a bus.
            busNames = cellstr([this.BusName]);
            busName = unique(busNames);
            % a valid spec should only have one unique bus name
            val = ~(numel(busName) > 1 || isempty(busName));
            if val
                busName = busName{1};
            else
                busName = '';
            end
        end
        
        function isCompat = isCompatible(obj,obj2)
            % check if data specs SPEC1 and SPEC2 have compatible size,
            % values, data type, limit, etc. Ignore differences in Name, Description
            % do not support matrix data specs
            
            if isa(obj2,'rl.util.RLDataSpec')
                isCompat = true;
                if isvector(obj) && isvector(obj2) && isequal(numel(obj),numel(obj2))
                    for ct = 1:numel(obj)
                        if ~(isa(obj(ct),class(obj2(ct))) && isCompatible_(obj(ct),obj2(ct)))
                            % return FALSE if individual elements are of
                            % different class or not compatible
                            isCompat = false;
                            break
                        end
                    end
                else
                    % not support matrix data specs
                    isCompat = false;
                end
            else
                error(message('rl:general:errInvalidCompareInputs'))
            end
        end
        
        function [BatchSize,SequenceLength] = inferDataDimension(obj, Data)
            % Return BatchSize and SequenceLength by comparing the data
            % dimension with the data spec. Convention:
            % DataSpecDim x BatchSize x SequenceLength
            % NOTE: for now only be called by representation class
            
            DataSz = size(Data);
            InputSz = obj.Dimension;
            NDimsDiff = numel(DataSz) - numel(InputSz);
            switch NDimsDiff
                case 2
                    BatchSize = DataSz(end-1);
                    SequenceLength = DataSz(end);
                case 1
                    BatchSize = DataSz(end);
                    SequenceLength = 1;
                case 0
                    BatchSize = 1;
                    SequenceLength = 1;
                otherwise
                    error(message('rl:agent:errRepIncorrectInputDim'));
            end
        end
    end
    
    %%
    methods (Access = public, Sealed)       
        function sample = usample(obj,BatchSize)
            % usample uniformly samples from the space
            if nargin == 1
                BatchSize = 1;
            elseif ~(isnumeric(BatchSize) && isscalar(BatchSize) && ...
                    isreal(BatchSize) && (BatchSize>0) && (rem(BatchSize,1) == 0))
                % Validate BatchSize
                error(message('rl:general:errValidateBatchSize'))
            end
            for ct = 1:length(obj)
                sample{1,ct} = usample_(obj(ct),BatchSize);
            end
        end
        
        function isValid = isDataValid(obj,Data)
            % isDataValid Validates that the data is consistent with data
            % specifications
            isValid = false;
            if ~iscell(Data)
                Data = {Data};
            end
            for ct = 1:numel(Data)
                isValid(ct) = isDataValid_(obj(ct),Data{ct});
            end
            isValid = all(isValid);
        end
    end
    
    %%
    methods (Abstract, Access = protected)
        % returns sample of the data spec
        sample = usample_(obj)
        % check if the input data is consistent with data specification
        b = isDataValid_(Data);
        % check if 2 scalar specs are consistent
        isCompat = isCompatible_(obj,obj2);
    end
end

