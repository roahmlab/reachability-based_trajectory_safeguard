classdef rlTableRepresentation < rl.util.rlAbstractRepresentation
    %rlTableRepresentation Tabular representation for Value and Q Tables
    %  Value Table V(s)  
    %       V = rlTable(ObservationInfo)
    %       obj = rlTableRepresentation(V)
    %  Q-Table Q(s,a)
    %       Q = rlTable(ObservationInfo,ActionInfo)
    %       obj = rlTableRepresentation(Q)
    %
    %  ObservationInfo and ActionInfo must be a scalar rlFiniteSetSpec
    %
    %  See also: rlFiniteSetSpec
    
    %   Copyright 2018-2019 The MathWorks, Inc.
    
    properties
        Table
    end
    
    methods
        function obj = rlTableRepresentation(Table,varargin)
            %rlTableRepresentation Construct an instance of the tablular
            %  represntation
            %  Value Table Representation
            %       obj = rlTableRepresentation(rlTable(ObservationInfo))
            %  Q Table Representation
            %       obj = rlTableRepresentation(rlTable(ObservationInfo,ActionInfo))            
            
            narginchk(1, 2)
            obj = obj@rl.util.rlAbstractRepresentation(Table,varargin{:});
            obj.Table = Table;
            [obj.ObservationInfo,obj.ActionInfo] = getInfo(Table);
            % REVISIT: do we need to set this.Validated = true?
            obj.Validated=true;
        end
        
        % Set the loss function
        function this = setLoss(this,CurrentLoss)
            CurrentLoss = validatestring(CurrentLoss,{'mse','cte'},'','Loss');
        end
        
        % Get computation model
        function mdl = getModel(this)
            mdl = this.Table;
        end
    end
    
    %% Implementation of Public Abstract Methods
    methods (Hidden)
        
        % Compute gradients of representation output wrt learnable parameters or inputs
        function[grad,parameter] = gradient(obj,dy,dx,Data,varargin)
            %% This function only supports the following pairs and API
            % Required: (dOutput/dParameter,dOutput/dAction,dOutput/dObservation,dLoss/dParameter)
            % dy is 'output' or 'loss'
            % dx is 'parameter','observation' or 'action'
            % inputValues is the values passed into the representation
            % [grad,parameter] = gradient(obj,'loss','parameters',inputValues,target)
            % [grad,parameter] = gradient(obj,'output','parameters',inputValues,'initialGradients',initgradients)% used for DDPG
            %  
            
            parameter = [];
            % REVISIT: should process data be on the abstract class
            %Data = cellfun(@(x) processData(obj,x),Data,'UniformOutput',false);
            
            if strcmpi(dy,'output')
                switch dx
                    case 'parameter'
                        grad = gradientWithRespectToParameters(obj.Table,Data);
                    case 'action'
                        error(message('rl:agent:errUnsupportedGradient','rlTable'))
                    case 'observation'
                        error(message('rl:agent:errUnsupportedGradient','rlTable'))
                    otherwise
                        error(message('rl:agent:errInvalidSyntaxGradient'))
                end
            elseif strcmpi(dy,'loss') && strcmpi(dx,'parameter')
                grad = gradientOfLossWithRespectToParameters(obj.Table,Data,varargin{1});
            else
                error(message('rl:agent:errInvalidSyntaxGradient'))
            end

        end
        
        % Step learnable parameters with gradient information or loss
        function obj = step(obj,gradient,~)
            %%
            Value = getLearnableParameterValues_(obj);
            Value{1} = Value{1} -0.5*obj.Options.LearnRate*gradient;
            obj = setLearnableParameterValues_(obj,Value);
            
        end
        
        % Fit learnable parameters with input and target data
        function obj = fit(obj,Data,Target)
            %% Data is the input to the network
            % Target is the target output of the network
            [gradVal] = gradient(obj,'loss','parameter',Data,Target);
            obj = step(obj,gradVal);
        end
        
        %% Copy representation
        function Representation2 = copy(obj)
            Representation2 = obj;
        end
    end
    
    %% Implementation of Protected Abstract Methods
    methods (Access=protected)
        function argStruct = generateEvaluateFunction_(this,argStruct)
            % generate the evaluate function for the layer representation
            outputstr = argStruct.EvaluateOutputString;
            inputstr  = argStruct.EvaluateInputString;
            policyname = argStruct.PolicyName;
            matfilename = argStruct.MATFileName;
            localfcnstr = argStruct.LocalFunctionString;
            evalfcnname = argStruct.EvaluateFunctionName;
            outputstr_y = rl.codegen.handleMultipleOutputStrings(outputstr);

            %% Extract vars and make sure the DAG can be constructed
            oinfo = this.ObservationInfo;
            ainfo = this.ActionInfo;
            
            if ~isscalar(ainfo) || ~isscalar(oinfo)
                error(message('rl:general:RepLayerCodeGenMIMONotSupported'));
            end
            
            %% Extract the DAG
            tbl = this.Table.Table; %#ok<NASGU>
            
            %% Save the network to disk
            eval(sprintf('%s = tbl;',policyname));
            save(matfilename,policyname);
            
            %% generate the local function
            load2persistentstr = sprintf([...
                'persistent %s\n',...
                'if isempty(%s)\n',...
                '\ts = coder.load(''%s'',''%s'');\n',...
                '\t%s = s.%s;\n',...
                'end'],...
                policyname,policyname,matfilename,policyname,policyname,policyname);
            [oelementstr,~,~] = rl.codegen.generateFiniteElementStrings(oinfo);
            [aelementstr,~,~] = rl.codegen.generateFiniteElementStrings(ainfo);
            localevalfcnstr = sprintf([...
                'function %s = %s(%s)\n',...
                '%s\n',...
                'actionSet = %s;\n',...
                'observationSet = %s;\n',...
                'actionIndex = rl.codegen.getElementIndex(actionSet,action);\n',...
                'observationIndex = rl.codegen.getElementIndex(observationSet,%s);\n',...
                '%s = %s(observationIndex,actionIndex);\n',...
                'end'],...
                outputstr_y,evalfcnname,inputstr,...
                load2persistentstr,...
                aelementstr,...
                oelementstr,...
                argStruct.InputString,...
                outputstr_y,policyname);
            
            %% attach it to the arg struct
            if isempty(localfcnstr)
                localfcnstr = localevalfcnstr;
            else
                localfcnstr = sprintf('%s\n%s',localfcnstr,localevalfcnstr);
            end
            argStruct.LocalFunctionString = localfcnstr;
        end
        function Value = evaluate_(obj,Data)
            % Get Value from Table
            Value = evaluate(obj.Table,Data);
        end
        
        function Value = getLearnableParameterValues_(obj)
            % Get Table Values
            Value = {obj.Table.Table};
        end
        
        function obj = setLearnableParameterValues_(obj,Value)
            % Set Table Values
            obj.Table.Table = Value{1};
        end
        
        function validateLearnableParameterSize(obj,Value)
            % No op as this will be handled at the rlTable level
        end
        
        function sz = getSize_(obj,identifier)
            switch identifier
                case 'observation'
                        sz{1,1} = size(getElementValue(obj.ObservationInfo,1));
                case 'action'
                    if isempty(obj.ActionInfo)
                        sz = {};
                    else
                        sz{1,1} = size(getElementValue(obj.ActionInfo,1));
                    end
                case 'output'
                    sz{1,1} = [1,1];
                case 'input'
                    % REVISIT
                    sz{1,1} = size(getElementValue(obj.ObservationInfo,1));
                    if ~isempty(obj.ActionInfo)
                        sz{1,2} = size(getElementValue(obj.ActionInfo,1));
                    end
                    
                case 'parameter'
                    sz = {size(obj.Table.Table)};
            end
        end
    end
end
