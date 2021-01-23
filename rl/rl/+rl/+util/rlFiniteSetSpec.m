classdef  rlFiniteSetSpec < rl.util.RLDataSpec
    % rlFiniteSetSpec: Creates an object to define the finite set of actions or observations.
    %
    %   spec = rlFiniteSetSpec(Elements) creates data spec with the discrete set
    %   defined by the elements Elements.
    %
    %   For example:
    %     spec = rlFiniteSetSpec([0;1;2;3])
    %     spec = rlFiniteSetSpec({[0,1];[1,1];[1,2];[1,3]})
    %
    %   Use "usample" to create a sample from this specification. For example:
    %      s = usample(spec)
    
    % Copyright 2018-2019 The MathWorks Inc.
    
    properties
        Elements 
    end
    
    properties (Access = protected,Dependent)
        NumElements
    end
    
    methods
        function obj = rlFiniteSetSpec(ElementList)
            % rlFiniteSetSpec creates a discrete data space
            %
            % obj = rl.util.rlFiniteSetSpec([0;1;2;3])
            % obj = rl.util.rlFiniteSetSpec({[0,1];[1,1];[1,2];[1,3]})
            
            narginchk(1,1)
            % ElementList must be a vecotr or a cell array with numeric
            % entries and same size

            b  = iscell(ElementList) && isnumeric(ElementList{1});
                    
            if ~(isvector(ElementList) && (isnumeric(ElementList)) || b)
                error(message('rl:general:errFiniteSetInvalidInput'))
            end
            
            % REVISIT specifying dimension size
            [elementSize,elementDataType] = getElementInfo(ElementList);
            obj = obj@rl.util.RLDataSpec(elementSize,elementDataType);
            
            % REVISIT error checking and if no Element list is set
            obj.Elements = ElementList(:);
           
        end
        
        function elementMat = getElementIndicationMatrix(obj,batchElement,batchSize)
            % return a T/F matrix with size (numElement,batchSize)
            % elementMat(i,j) = true indicates element jth of batchElement matches
            % element ith of the data spec
            % For example, if the data spec has 2 elements (1 and 2) and
            % batchElement has batchSize of 3 (2,1,2), 
            % then elementMat = [0,1,0;
            %                    1,0,1]
            % batchElement is always a single-datatype matrix (output of DLT neural networks) 
        
            num = obj.NumElements; % number of elements
    
            % Initialize element indication matrix
            elementMat = false(num,batchSize); 
            % Get element size 
            [elemVectorLength, elemNonSingletonDim] = max(obj.Dimension);
            % Get element values (in single)
            elements = obj.Elements;
            % Force data spec's element list to be column vector
            elements = reshape(elements,1,[]);
            % Concat all element in the element list to a matrix with dimension = (elementSize,batchSize)
            if iscell(elements) && isnumeric(elements{1})
                % Force each element in the data spec's element list to be column vector
                if elemNonSingletonDim == 2
                    elements = cellfun(@transpose,elements,'UniformOutput',false);
                end
                allElementValues = cell2mat(elements);
            else
                % Each element in the data spec's element list is scalar
                % (list is column vector)
                allElementValues = elements;
            end

            allElementValues = reshape(allElementValues,elemVectorLength,[]);
            elementValuesFromSpec = single(allElementValues)'; % single is required by DLT            
            elementBatchFromInput = reshape(batchElement{1},elemVectorLength,[])';

            % Each row in elementValuesFromSpec contains action.
            % elementMat contains one-hot-vector in each row that indicates which
            % action is selected. The following part loops over actions
            % (num) and finds the index of selected actions.
            for ct = 1:num
                elementMat(ct,:) = reshape(all(elementBatchFromInput==elementValuesFromSpec(ct,:),2),size(elementMat(ct,:)));
            end
        end
        
        %% Get methods
        function num = get.NumElements(obj)
            % Returns the current number of elements 
            num = numel(obj.Elements);
        end
        
        function num = getNumberOfElements(obj)
            % Number of Elements
            num = [obj.NumElements];
            num = reshape(num,size(obj));
        end
        
        function Value = getElementValue(obj,Idx)
            % Value of element idx
            % numel(obj) == numel(Idx)
            
            NumChannel = numel(obj);
            if NumChannel < 2
                ElementList = obj.Elements;
                if iscell(ElementList)
                    Value = ElementList{Idx};
                else
                    Value = ElementList(Idx);
                end
            else
                Value = cell(1,NumChannel);
                for ct = 1:NumChannel
                    ElementList = obj(ct).Elements;
                    if iscell(ElementList)
                        Value{ct} = [ElementList{Idx(ct,:)}];
                    else
                        Value{ct} = ElementList(Idx(ct,:));
                    end
                end
            end
        end
        
        function idx = getElementIndex(obj,Data)
            ElementList = obj.Elements;
            nData = size(Data,4);
            idx = zeros(nData,1);
            for ct = 1:nData
                D = Data(:,:,:,ct);
                if iscell(ElementList)
                    cidx = [];
                    for ct2 = 1:numel(ElementList)
                        if isequal(D,ElementList{ct2})
                            cidx = ct2;
                            break
                        end
                    end
                else                   
                    cidx = find(D == ElementList);
                end
                if isempty(cidx)
                    error(message('rl:general:errFiniteSetInvalidData'))
                else
                    idx(ct,1) = cidx;
                end
            end
        end
        
        %% Set methods
        function obj = set.Elements(obj,ElementList)
            obj.Elements = ElementList;
            if ~isempty(obj.Elements)
                [obj.Dimension,obj.DataType] = getElementInfo(ElementList);
            end
        end
        
    end
    
    
    %% Implementation of absract methods
    methods (Access = protected)
        function samples = usample_(obj,BatchSize)
            % usample uniformly samples from the space
            samples = [];
            for ct = 1:BatchSize
                % Revisit cat dimension
                samples = cat(2,samples,getElementValue(obj,randi(obj.NumElements)));
            end
        end
        
        function isValid = isDataValid_(obj,Data)
            % isDataValid Validates that the data is consistent with data
            % specifications (data type, dimension, range)
            isValid = isa(Data,obj.DataType) && isequal(size(Data), obj.Dimension);
            if iscell(obj.Elements) && ~iscell(Data)
                isValid = isValid && any(cellfun(@(x) isequal(Data,x),obj.Elements));
            else
                isValid = isValid && ismember(Data,obj.Elements);
            end
        end
        
        function isCompat = isCompatible_(obj,obj2)
            % isCompatible_ checks if 2 scalar finite set specs are consistent
            
            isCompat = isequal(obj.Elements,obj2.Elements) && ...
                isequal(obj.Dimension,obj2.Dimension) && ...
                isequal(obj.DataType,obj2.DataType);
        end
        
    end
end

%% Local function: get size and dataType of element
function [elementSize,elementDataType] = getElementInfo(ElementList)
if iscell(ElementList)
    firstElement = ElementList{1};
    if numel(ElementList) > 1 && ~isa(firstElement,'char')
        % check all input element have the same dimension
        if ~all(cellfun(@(element) isequal(size(firstElement), size(element)) , ElementList(2:end)))
            error(message('rl:general:errFiniteSetElementMismatchDim'));
        end
        % check all input element have the same data type
        if ~all(cellfun(@(element) isa(element, class(firstElement)) , ElementList(2:end)))
            error(message('rl:general:errFiniteSetElementMismatchType'));
        end
    end
else
    firstElement = ElementList(1);
end
elementSize = size(firstElement);
elementDataType = class(firstElement);
end