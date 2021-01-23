classdef GaussianActionNoise < rl.option.AbstractActionNoise
    % GAUSSIANACTIONNOISE: Creates options for noise model used 
    % in deterministic policy gradient agent. You can also create 
    % GaussianActionNoise objects for your own custom agents.
    %
    % Supported Options are:
    %    Mean: mean of noise model
    %    Variance: variance of noise model
    %    VarianceDecayRate: decay rate of variance
    %    VarianceMin: Minimum value of variance
    %    LowerLimit: lower limit of noise samples
    %    UpperLimit: upper limit of noise samples
    % 
    % Noise is sampled according the following formula:
    % 
    % x = Mean + rand(ActionSize).*Variance
    % x = min(max(x,LowerLimit),UpperLimit);
    % 
    % Examples:
    % obj = rl.option.GaussianActionNoise;
    % obj = rl.option.GaussianActionNoise('Mean',0.5);
    % obj = rl.option.GaussianActionNoise('Mean',0.5,...
    %                                     'Variance',0.4);
    % obj = rl.option.GaussianActionNoise('Variance',0.1,...
    %                                     'Clip',0.5);
    %
    % An GaussianActionNoise object is automatically created within 
    % an rlTD3AgentOptions object. In this case, specify noise model options
    % using dot notation.
    %
    % opt = rlTD3AgentOptions;
    % opt.ExplorationNoiseOptions.Variance = 0.2;
    
    % Copyright 2019 The MathWorks Inc.
    
    properties       
        % Noise mean
        Mean
        % Decay rate of the variance (0 -> no decay, 1-> immediate decay)
        VarianceDecayRate
    end
    
    properties (Dependent)        
        % Variance of the noise model
        Variance
        % Minimum value of variance
        VarianceMin
        % Lower bound of noise samples
        LowerLimit
        % Upper bound of noise samples
        UpperLimit
    end
    
    properties (Access = private)       
        % Variance of the noise model
        Variance_
        % Minimum value of variance
        VarianceMin_
        % Lower bound of noise samples
        LowerLimit_
        % Upper bound of noise samples
        UpperLimit_
    end
    
    methods
        function obj = GaussianActionNoise(varargin)
            parser = inputParser;
            addParameter(parser,'Mean',0)
            addParameter(parser,'Variance',0.2)
            addParameter(parser,'VarianceDecayRate',0)
            addParameter(parser,'VarianceMin',0.01)
            addParameter(parser,'LowerLimit',-Inf)
            addParameter(parser,'UpperLimit',Inf)
            parse(parser,varargin{:});

            obj.Mean = parser.Results.Mean;
            obj.Variance = parser.Results.Variance;
            obj.VarianceDecayRate = parser.Results.VarianceDecayRate;
            obj.VarianceMin = parser.Results.VarianceMin;
            obj.LowerLimit = parser.Results.LowerLimit;
            obj.UpperLimit = parser.Results.UpperLimit;
        end
        
        function obj = set.Mean(obj,Value)
            PropertyName = 'Mean';
            ValFcn = @(Value) validateattributes(Value,{'numeric'},{'real','nonempty','finite','nonnan'},'',PropertyName);
            checkParam(obj,PropertyName,ValFcn,Value);
            obj.Mean = Value;
        end
        
        function obj = set.VarianceDecayRate(obj,Value)
            PropertyName = 'VarianceDecayRate';
            ValFcn = @(Value) validateattributes(Value,{'numeric'},{'real','nonnegative','<=',1,'nonempty','nonnan'},'',PropertyName);
            checkParam(obj,PropertyName,ValFcn,Value);
            obj.VarianceDecayRate = Value;
        end

        function obj = set.Variance(obj,Value)
            PropertyName = 'Variance';
            ValFcn = @(Value) validateattributes(Value,{'numeric','cell'},{'real','nonnegative','nonempty','finite','nonnan'},'',PropertyName);
            checkParam(obj,PropertyName,ValFcn,Value);
            if ~isempty(obj.Variance) || ~isempty(obj.VarianceMin)
                checkParameterValues(obj,Value,'<','VarianceMin','warning','rl:agent:warnNoiseVarianceLtVarianceMin');
            end
            obj.Variance_ = Value;
        end
        
        function obj = set.VarianceMin(obj,Value)
            PropertyName = 'VarianceMin';
            ValFcn = @(Value) validateattributes(Value,{'numeric','cell'},{'real','nonnegative','nonempty','finite','nonnan'},'',PropertyName);
            checkParam(obj,PropertyName,ValFcn,Value);
            if ~isempty(obj.Variance) || ~isempty(obj.VarianceMin)
                checkParameterValues(obj,Value,'>','Variance','warning','rl:agent:warnNoiseVarianceLtVarianceMin');
            end
            obj.VarianceMin_ = Value;
        end

        function obj = set.LowerLimit(obj,Value)
            PropertyName = 'LowerLimit';
            ValFcn = @(Value) validateattributes(Value,{'numeric','cell'},{'real','nonempty','nonnan'},'',PropertyName);
            checkParam(obj,PropertyName,ValFcn,Value);
            if ~isempty(obj.LowerLimit) || ~isempty(obj.UpperLimit)
                checkParameterValues(obj,Value,'>','UpperLimit','error','rl:agent:errGaussianActionNoiseInvalidBounds');
            end
            obj.LowerLimit_ = Value;
        end
        
        function obj = set.UpperLimit(obj,Value)
            PropertyName = 'UpperLimit';
            ValFcn = @(Value) validateattributes(Value,{'numeric','cell'},{'real','nonempty','nonnan'},'',PropertyName);
            checkParam(obj,PropertyName,ValFcn,Value);
            if ~isempty(obj.LowerLimit) || ~isempty(obj.UpperLimit)
                checkParameterValues(obj,Value,'<','LowerLimit','error','rl:agent:errGaussianActionNoiseInvalidBounds');
            end
            obj.UpperLimit_ = Value;
        end
        
        function Options = get.Variance(this)
            Options = this.Variance_;
        end
        
        function Options = get.VarianceMin(this)
            Options = this.VarianceMin_;
        end
        
        function Options = get.LowerLimit(this)
            Options = this.LowerLimit_;
        end
        
        function Options = get.UpperLimit(this)
            Options = this.UpperLimit_;
        end
    end
    
    methods (Access = protected)
        function noiseProperties = getNoiseProperties(~)
            noiseProperties = {...
                    'Mean',...
                    'Variance',...
                    'VarianceDecayRate',...
                    'VarianceMin',...
                    'LowerLimit',...
                    'UpperLimit',...
                    };
        end
    end
end