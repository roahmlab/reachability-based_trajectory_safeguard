classdef rlLinearBasisRepresentation < rl.util.rlModelRepresentation
    %rlLinearBasisRepresentation Linear Basis representation for Value and
    %  Q functions and policies.
    %
    %  Value Function V(s)  
    %       obj = rlLinearBasisRepresentation(myValueFcn,w0,OINFO)
    %
    %  Q-Table Q(s,a)
    %       obj = rlLinearBasisRepresentation(myQFcn,w0,{OINFO,AINFO})
    %
    %  Policy c(s)
    %       obj = rlLinearBasisRepresentation(myPolicyFcn,w0,OINFO,AINFO)
 
    %   Copyright 2018-2019 The MathWorks, Inc.
    
    properties
        % Function handle to basis function
        FcnHandle
    end
    
    properties (Access = private)       
        % Representation type V, Q, P
        Type
    end
    
    methods
        function obj = rlLinearBasisRepresentation(BasisFcn,w0,varargin)
            %rlLinearBasisRepresentation Construct an instance of the
            %  linear basis representation represntation
            %
            %  Value Table Representation
            %       obj = rlLinearBasisRepresentation(ObservationInfo)
            %  Q Table Representation
            %       obj = rlLinearBasisRepresentation(ObservationInfo,ActionInfo) 
            
            ix = find(cellfun(@(x) isa(x,'rl.option.rlRepresentationOptions'),varargin),1);
            if isempty(ix)
                Options = {};              
            else
                Options = varargin(ix);
                varargin(ix) = [];
            end        
            [model,args,Type] = localLinearBasis(BasisFcn,w0,varargin{:});
            args = [args,Options];
            obj = obj@rl.util.rlModelRepresentation(model,args{:});     
            obj.FcnHandle = BasisFcn;
            obj.Type = Type;
            
            
        end
        
        function s = saveobj(obj)
            s.FcnHandle = obj.FcnHandle;
            s.Params = obj.getLearnableParameterValues;
            s.ObsInfo = obj.ObservationInfo;
            s.ActInfo = obj.ActionInfo;
            s.Type = obj.Type;
            s.Options = obj.Options;
            Loss = obj.adModel.getLoss;
            s.LossName = Loss.Name;

        end
        
        % Get computation model
        function mdl = getModel(this)
            % REVISIT: Consider return a function handle with learned
            % weights. Might reuse generateEvaluateFunction_
            % Currently output a struct of function handle and learnable 
            % parameters for simplicity.
            mdl = {this.FcnHandle getLearnableParameterValues(this)};
        end
       
    end 
    
    methods (Hidden)
       
        %% Copy representation
        function Representation2 = copy(Representation1)
            % Representation is value object
            Representation2 = Representation1.loadobj(saveobj(Representation1));
        end
        
    end
    
    methods (Static)
        function obj = loadobj(s)
            % Recontstruct object
            switch s.Type
                case "V"
                    obj = rlRepresentation(s.FcnHandle, s.Params{1}, s.ObsInfo,s.Options);
                case "Q"
                    obj = rlRepresentation(s.FcnHandle, s.Params{1}, {s.ObsInfo,s.ActInfo},s.Options);
                case "P"
                    obj = rlRepresentation(s.FcnHandle, s.Params{1}, s.ObsInfo, s.ActInfo,s.Options);
            end
            obj = obj.setLoss(s.LossName);
        end
    end
    
    methods (Access=protected)
        function argStruct = generateEvaluateFunction_(this,argStruct)
            % generate the evaluate function for the layer representation
            outputstr = argStruct.EvaluateOutputString;
            inputstr  = argStruct.EvaluateInputString;
            %             policyname = argStruct.PolicyName;
            %             matfilename = argStruct.MATFileName;
            localfcnstr = argStruct.LocalFunctionString;
            evalfcnname = argStruct.EvaluateFunctionName;
            outputstr_y = rl.codegen.handleMultipleOutputStrings(outputstr);
            
            %%
            fcn = this.FcnHandle;
            wcell = getLearnableParameterValues(this);
            w = wcell{1};
            
            fcnstr = func2str(fcn);
            wstr = mat2str(w);
            
            str1 = sprintf([...
                'fcn = %s;\n',...
                'w = %s;\n'],...
                fcnstr,...
                wstr);
            % REVISIT: v1 assume single loss
            modelOutputs = getOutputs(this.adModel,1);
            if isa(modelOutputs,'rl.internal.ad.ops.Softmax') || isa(modelOutputs,'rl.ad.ops.Softmax')
                % add softmax to the generated code
                str2 = sprintf('%s = rl.codegen.softmax(w''*fcn(%s));\n',...
                    outputstr_y,inputstr);
            else
                str2 = sprintf('%s = w''*fcn(%s);\n',...
                    outputstr_y,inputstr);
            end
            localevalfcnstr = sprintf([...
                'function %s = %s(%s)\n',...
                '%s%s',...
                'end'],...
                outputstr_y,evalfcnname,inputstr,...
                str1,str2);
            
            %% attach it to the arg struct
            if isempty(localfcnstr)
                localfcnstr = localevalfcnstr;
            else
                localfcnstr = sprintf('%s\n%s',localfcnstr,localevalfcnstr);
            end
            argStruct.LocalFunctionString = localfcnstr;
        end
        
        function this = updateOptions(this,opt)
            % convert representation option to local options
            this.LocalOptions = rl.util.rlModelRepresentation.createOptimizerOptions(opt,1);
            if this.Validated
                setOptimizer(this.adModel,this.LocalOptions);
            end
        end
        
    end
            
    
end

%% Local Functions
function [model,args,Type] = localLinearBasis(BasisFcn,w0,InputInfo,ActInfo)
% rlLinearBasis Create a linear basis for a critic or actor
%
%  V = rlLinearBasis = rlLinearBasis(BasisFcn, w0, ObsInfo)
%      creates f(s) = w'*B(s) Value function
%
%  Q = rlLinearBasis = rlLinearBasis(BasisFcn, w0, {ObsInfo, ActInfo})
%      creates f(s,a) = w'*B(s,a) Q function
%
%  P = rlLinearBasis = rlLinearBasis(BasisFcn, w0, ObsInfo, ActInfo)
%      creates f(s,a) = w'*B(s) Policy function
%
%  where
%       BasisFcn is a function handle that returns a column vectorn B for the basis
%       w0 is the inital paramters
%       ObsInfo is the observation information for the agent
%       ActInfo is the action information for the agent
narginchk(3,4)

 if nargin == 3 
     % Type = Q or V
     if iscell(InputInfo) 
         if numel(InputInfo) == 1
             ObsInfo = InputInfo{1};
             ActInfo = [];
             Type = "V";
         elseif numel(InputInfo) == 2
             ObsInfo = InputInfo{1};
             ActInfo = InputInfo{2};
             validateattributes(ActInfo,{'rl.util.RLDataSpec'},{'vector'},'','ActInfo')
             Type = "Q";
         end
     else
         ObsInfo = InputInfo;
         ActInfo = [];
         Type = "V";
     end
 end
 
 if nargin == 4
     % Type = P
     Type = "P";
     ObsInfo = InputInfo;
 end
 
 validateBasisFcn(BasisFcn,w0,Type,ObsInfo,ActInfo)
 
 validateattributes(ObsInfo,{'rl.util.RLDataSpec'},{'vector'},'','ObsInfo')
     
[s,ObsInfo] = localInfo2adInput(ObsInfo,"O");
w = rl.internal.ad.Var('Name', 'w', 'Value', w0);

switch Type
    case 'V'
        out = w'*feval(BasisFcn,s{:});
        model = rl.internal.ad.Model({s{1,:}}, out);
        args = {ObsInfo,'Observation',cellstr([ObsInfo.Name])};
    case 'Q'
        [a,ActInfo] = localInfo2adInput(ActInfo,"A");
        out = w'*feval(BasisFcn,s{:},a{:});
        model = rl.internal.ad.Model({s{1,:},a{1,:}}, out);
        args = {ObsInfo,'Observation',cellstr([ObsInfo.Name]),ActInfo,'Action',cellstr([ActInfo.Name])};
    case 'P'
        [~,ActInfo] = localInfo2adInput(ActInfo,"A");
        out = w'*feval(BasisFcn,s{:});
        out.Name = ActInfo.Name;
        model = rl.internal.ad.Model({s{1,:}}, out);
        args = {ObsInfo,'Observation',cellstr([ObsInfo.Name]),ActInfo,'Action',cellstr([ActInfo.Name])};
end

end



function [adin,info] = localInfo2adInput(info,Name)
   % Convert info vector into input variables
   for ct = 1:numel(info)
       value = usample(info(ct));
       if isempty(info(ct).Name)
           info(ct).Name = strcat(Name,num2str(ct));
       end
       adin{1,ct} = rl.internal.ad.Input(size(value{1}),'Value',value{1},'Name',char(info(ct).Name),'Precision','double');
   end

end


function validateBasisFcn(BasisFcn,w0,Type,ObsInfo,ActInfo)
switch Type
    case 'V'
        ObsData = usample(ObsInfo);
        try
            B = feval(BasisFcn,ObsData{:});
        catch
            error('The BasisFcn input signature does not match.')
        end
        try
            F = w0'*B;
        catch
            error('The weights w are not compatible with the basis function.')
        end
        if ~isscalar(F)
            error('The output of the value function V=w''*B must be a scalar.')
        end
        
        
    case 'Q'
        ObsData = usample(ObsInfo);
        ActData = usample(ActInfo);
        try
            B = feval(BasisFcn,ObsData{:},ActData{:});
        catch
            error('The BasisFcn input signature does not match.')
        end
        try
            F = w0'*B;
        catch
            error('The weights w are not compatible with the basis function.')
        end
        if ~isscalar(F)
            error('The output of the Q function Q=w''*B must be a scalar.')
        end
        
    case 'P'
        ObsData = usample(ObsInfo);
        try
            B = feval(BasisFcn,ObsData{:});
        catch
            error('The BasisFcn input signature does not match.')
        end
        try
            F = w0'*B;
        catch
            error('The weights w are not compatible with the basis function')
        end
        
end

end


