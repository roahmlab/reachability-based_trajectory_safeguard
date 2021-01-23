classdef  rlNumericSpec < rl.util.RLDataSpec
    % rlNumericSpec N-dimensional data class for action and observation space
    % definitions 
    
    % Copyright 2018-2019 The MathWorks Inc.
    
    % other names:
    %   RLDataBounded
    %   RLDataReal
    %   RLDataNumeric
    
    properties
        % Lower limit for the data space. Must be a scalar value or a
        % matrix of same size of the data space. If value is scalar it will
        % be applied to all entries in the data space.
        LowerLimit
        
        % Upper limit for the data space. Must be a scalar value or a
        % matrix of same size of the data space. If value is scalar it will
        % be applied to all entries in the data space.
        UpperLimit 
    end
    
    properties (Access = private)
        % Version indicator, increase this value by 1 when modyfing loadobj
        Version = 1
    end
    
    methods
        function obj = rlNumericSpec(Dimension, varargin)
            % RLDataTensor creates a N-Dimensional data space
            %
            %  obj = rl.util.rlNumericSpec([160,210,3], ...
            %       'LowerLimit', 0, 'UpperLimit', 1, 'DataType', 'single');
            
            validateattributes(Dimension, {'numeric'}, {'nonsparse', 'row', 'finite'}, mfilename, 'Dimension', 1);
            
            % remove trailing 1 for dimension. 
            % e.g. [2 1 1 1] will turn to [2 1]
            Dimension = size(zeros(Dimension));
            
            % Parse varargin
            parser = inputParser;
            addParameter(parser,'LowerLimit',-inf);
            addParameter(parser,'UpperLimit',inf);
            addParameter(parser,'DataType',"double");
            parse(parser,varargin{:})
            
            % Call parent contstructor
            obj = obj@rl.util.RLDataSpec(Dimension,parser.Results.DataType);
            
            % Get Lower Limit and Upperlimit
            LL = parser.Results.LowerLimit;
            UL = parser.Results.UpperLimit;
            
            % Validate lower and upper limit
            % Limit sizes must be scalar or size of Dimension
            % Lower Limit must be less than Upper Limit
            if ~(isscalar(LL)|| isequal(size(LL),Dimension))
                error(message('rl:general:errUpperLowerBoundSize','LowerBound',mat2str(Dimension)))
            end
            
            if ~(isscalar(UL)|| isequal(size(UL),Dimension))
                error(message('rl:general:errUpperLowerBoundSize','UpperBound',mat2str(Dimension)))
            end
            
            if all(LL<=UL,'all')
                % Set upper and lower limits
                obj.LowerLimit = LL;
                obj.UpperLimit = UL;
            else
                error(message('rl:general:errUpperBoundGreaterThanLowerBouund'))
            end
            
            % Version indicator, increase this value by 1 when modyfing loadobj
            obj.Version = 2;
        end
    end
    
    methods (Hidden)
        function Data = saturate(obj, Data)
            % SATURATE saturates the data based on rlNumericData spec 
            % LowerLimit and UpperLimit. 
            %
            %  DATA = saturate(NUMERICSPEC, DATA)
            %   - If NUMERICSPEC is a scalar, DATA can be numeric or a cell
            %   - If NUMERICSPEC is a vector, DATA should be a cell array
            %   with the same size as NUMERICSPEC.
            %   - Extra elements in DATA will be truncated.
            
            try
                if iscell(Data)
                    for i = 1:numel(obj)
                        Data{i} = min(max(Data{i}, obj(i).LowerLimit), obj(i).UpperLimit);
                    end
                else
                    Data = min(max(Data, obj.LowerLimit), obj.UpperLimit);
                end
                Data = Data(1:numel(obj));  % truncate extra elements
            catch
                validateattributes(Data, {'numeric','cell'}, {'nonempty'}, '', 'Data');
                if iscell(Data)
                    validateattributes(Data, {'cell'}, {'vector','nonempty'}, '', 'Data');
                    if (numel(Data) ~= numel(obj))
                        error(message('rl:general:errNumericSaturateMismatchNumel'))
                    end
                else
                    if numel(obj) ~= 1
                        error(message('rl:general:errNumericSaturateMismatchNumel'))
                    end
                end
            end
        end
        
        function [Scale, Bias] = getTanhScaleBias(obj)
            % GETTANHSCALEBIAS gets the scale and bias to move the outputs 
            % of tanh() (bound [-1 1]) to the data spec bound of interest.
            %
            %  [SCALE, BIAS] = getTanhScaleBias(NUMERICSPEC)
            %   - If NUMERICSPEC is a scalar, SCALE and BIAS are 1x1 cells
            %   - If NUMERICSPEC is a vector, SCALE and BIAS are cell
            %   arrays with the same size as NUMERICSPEC.
            %   - Each element of SCALE/BIAS cell has the same size with
            %   NUMERICSPEC's bound size.
            
            NumNumericSpec = numel(obj);
            Scale = cell(NumNumericSpec,1);
            Bias  = cell(NumNumericSpec,1);
            
            for ct = 1:NumNumericSpec
                % UpperLimit and LowerLimit have same dimension or be scalar
                CurrentScale = (obj(ct).UpperLimit - obj(ct).LowerLimit) / 2;
                CurrentBias  = obj(ct).UpperLimit - CurrentScale;
                
                % no-op on element with Scale or Bias of Inf (default bound)
                CurrentScale(isinf(CurrentScale)) = 1;
                CurrentBias(isinf(CurrentBias)) = 0;
                
                Scale{ct} = CurrentScale;
                Bias{ct} = CurrentBias;
            end
        end
    end
    
    %% Implement Abstract Methods
    methods (Access = protected)      
        function samples = usample_(obj,BatchSize)
            % usample uniformly samples from the space
            
            if numel(obj.LowerLimit) == 1
                LL = obj.LowerLimit*ones(obj.Dimension);
            else
                LL = obj.LowerLimit;
            end

            if numel(obj.UpperLimit) == 1
                UL = obj.UpperLimit*ones(obj.Dimension);
            else
                UL = obj.UpperLimit;
            end
            
            % REVISIT: Handle infs
            LL(LL==-inf) = -10000;
            UL(UL==inf) = 10000;
            
            
            samples = [];
            for ct = 1:BatchSize
               samples = cat(numel(obj.Dimension)+1,samples,(UL-LL).*rand(obj.Dimension) + LL);
            end
            
        end
        
        function isValid = isDataValid_(obj,Data)
            % isDataValid Validates that the data is consistent with data
            % specifications (data type, dimension, range)
            isValid = isa(Data,obj.DataType) && isequal(size(Data), obj.Dimension) && ...
                all(obj.LowerLimit <= Data & Data <= obj.UpperLimit,'all');
        end
        
        function isCompat = isCompatible_(obj,obj2)
            % isSpecCompatible_ checks if 2 scalar numeric specs are consistent
            
            isCompat = isequal(obj.LowerLimit,obj2.LowerLimit) && ...
                isequal(obj.UpperLimit,obj2.UpperLimit) && ...
                isequal(obj.Dimension,obj2.Dimension) && ...
                isequal(obj.DataType,obj2.DataType);
        end
    end

    methods (Static)
        function obj = loadobj(s)
            if s.Version == 1
                % version 1 does not have Version property so will use
                % default value
                % In version 2,
                %   - remove trailing 1 for dimension
                obj             = rlNumericSpec(s.Dimension);
                obj.LowerLimit  = s.LowerLimit;
                obj.UpperLimit  = s.UpperLimit;
                obj.Name        = s.Name;
                obj.Description = s.Description;
                obj.DataType    = s.DataType;
                obj.BusName     = s.BusName;
            else
                obj = s;
            end
        end
    end
end