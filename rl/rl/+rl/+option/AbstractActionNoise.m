classdef AbstractActionNoise
    
    methods (Access = protected)
        function checkParam(this,name,valfcn,Value)
            % Util method to check data dimension of each property when set
            % Each property can be a scalar (expansion) or a cell array.
            % All properties have same number of cells (action channels)
            % Each cell (Each action channel) all properties must have same
            % number of elements and dimensions.
            
            if iscell(Value)
                cellfun(valfcn,Value);
            else
                valfcn(Value);
            end
            checkParameterConsistency(this,name,Value);
        end
        function checkParameterConsistency(this,name,val)
            % the parameters below must either have the same dimension or
            % be scalar
            if ~isscalar(val) || iscell(val)
                params_ = getNoiseProperties(this);
                params_ = setdiff(params_,{name});
                
                % make sure all the cells have the same number of elements
                % and dimensions
                cells_ = cellfun(@(x) iscell(this.(x)),params_);
                if any(cells_)
                    cparams_ = params_(cells_);
                    
                    % make sure the number of cells is consistent
                    if iscell(val)
                        nc = cellfun(@(x) numel(this.(x)),cparams_);
                        [nc,maxidx] = max(nc);
                        nv = numel(val);
                        if nv ~= nc
                            error(message('rl:agent:errOUNoiseInconsistentNumActions',name,mat2str(nv),cparams_{maxidx},mat2str(nc)));
                        end
                    end
                    
                    % make sure the dims for each element is consistent
                    fcn = @(v) cellfun(@(x)size(squeeze(x)),v,'UniformOutput',false);
                    if iscell(val)
                        v1 = fcn(val);
                    else
                        v1 = fcn({val});
                    end
                    for i = 1:numel(cparams_)
                        p = cparams_{i};
                        v2 = fcn(this.(p));
                        nv2 = numel(v2);
                        if numel(v1) ~= nv2
                            v1 = repmat(v1,1,nv2);
                        end
                        invalidIdx = find(cellfun(@(x,y) ~isequal(x,y),v1,v2),1);
                        if ~isempty(invalidIdx)
                            error(message('rl:agent:errOUNoiseInconsistentActionSizes',invalidIdx,name,mat2str(v1{invalidIdx}),invalidIdx,p,mat2str(v2{invalidIdx})));
                        end
                    end
                    params_ = params_(~cells_);
                end
                
                sizes_ = cellfun(@(x) size(squeeze(this.(x))),params_,'UniformOutput',false);
                numels_ = cellfun(@prod,sizes_);
                [~,maxidx] = max(numels_);

                maxsize = sizes_{maxidx};
                valsize = size(squeeze(val));
                if ~isequal(valsize,maxsize) && prod(maxsize) > 1
                    error(message('rl:agent:errOUNoiseInconsistentParamSizes',name,mat2str(valsize),params_{maxidx},mat2str(maxsize)));
                end
            end
        end
        function checkParameterValues(obj,Value,Operator,PropName,Type,errID)
        % CHECKPARAMETERVALUES Throw error/warning errID if the condition
        % Value <Operator> PropName is true.
        %
        % Value : Value to compare
        % Operator : logic
        % PropName: Property to compare against
        % Type : error or warning
        % errID : error ID
        %
            Prop = obj.(PropName);  % Get the property by its name
            
            if isnumeric(Value) && isnumeric(Prop)
                % convert to vectors for comparison
                ValueVector = reshape(Value,1,numel(Value));
                PropVector = reshape(Prop,1,numel(Prop));
                checkValues();
            elseif isnumeric(Value) && iscell(Prop)
                for i = 1:numel(Prop)
                    % convert to vectors for comparison
                    PropVector = reshape(Prop{i},1,numel(Prop{i}));
                    if isscalar(Value)   % Support implicit expansion
                        ValueVector = Value;
                    else
                        ValueVector = reshape(Value{i},1,numel(Value{i}));
                    end
                    checkValues();
                end
            elseif iscell(Value)
                for i=1:numel(Value)
                    % convert to vectors for comparison
                    ValueVector = reshape(Value{i},1,numel(Value{i}));
                    if isscalar(Prop)   % Support implicit expansion
                        PropVector = Prop;
                    else
                        PropVector = reshape(Prop{i},1,numel(Prop{i}));
                    end
                    checkValues();
                end
            end
            
            % Nested function for reuse
            function checkValues()
                switch Operator
                    case '<'
                        if any(ValueVector < PropVector)
                            if strcmp(Type,'warning')
                                warning(message(errID));
                            elseif strcmp(Type,'error')
                                error(message(errID));
                            end
                        end
                    case '>'
                        if any(ValueVector > PropVector)
                            if strcmp(Type,'warning')
                                warning(message(errID));
                            elseif strcmp(Type,'error')
                                error(message(errID));
                            end
                        end
                end
            end
            
        end
    end
    methods (Access = protected)
        % Return all the properties of the noise to check for consistency
        noiseProperties = getNoiseProperties(this);
    end
end