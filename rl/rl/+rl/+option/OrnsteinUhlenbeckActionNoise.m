classdef OrnsteinUhlenbeckActionNoise < rl.option.AbstractActionNoise
    % ORNSTEINUHLENBECKACTIONNOISE: Creates options for noise model used 
    % in DDPG agent. You can also create OrnsteinUhlenbeckActionNoise
    % objects for your own custom agents.
    %
    % Supported Options are:
    %    InitialAction: initial value of action for noise model
    %    Mean: mean of noise model
    %    MeanAttractionConstant: constant that decides how quickly output 
    %    of noise model is attracted to Mean
    %    Variance: variance of noise model
    %    VarianceDecayRate: decay rate of variance
    %    SampleTime: sample time for noise model update
    % 
    % Noise model is updated according the following formula:
    % 
    % x(k) = x(k-1) + MeanAttractionConstant.*(Mean - x(k-1)).*SampleTime...
    %        + Variance.*randn(size(Mean)).*sqrt(SampleTime)
    % 
    % Examples:
    % obj = rl.option.OrnsteinUhlenbeckActionNoise;
    % obj = rl.option.OrnsteinUhlenbeckActionNoise('Mean',0.5);
    % obj = rl.option.OrnsteinUhlenbeckActionNoise('Mean',0.5,...
    %                                               'Variance',0.4);
    % obj = rl.option.OrnsteinUhlenbeckActionNoise('Variance',0.4,...
    %                                              'VarianceDecayRate',1e-5);
    %
    % An OrnsteinUhlenbeckActionNoise object is automatically created within 
    % an rlDDPGAgentOptions object. In this case, specify noise model options
    % using dot notation.
    %
    % opt = rlDDPGAgentOptions;
    % opt.NoiseOptions.Mean = 0.5;
    
    % Copyright 2019 The MathWorks Inc.
    
    properties
        % Initial noise
        InitialAction
        % Noise mean
        Mean
        % Attraction parameter of the model to get back to the mean
        MeanAttractionConstant
        % Decay rate of the variance (0 -> no decay, 1-> immediate decay)
        VarianceDecayRate
    end
    properties (Dependent)
        % Variance of the noise model
        Variance
        % Minimum variance of the noise model
        VarianceMin
    end
    properties (Access = private)
        % Variance of the noise model
        Variance_
        
        % Minimum variance of the noise model
        VarianceMin_ = 0
    end
    properties (Hidden)
        % Sample time of the process, use -1 for inherited from Agent
        SampleTime
    end
    methods
        function obj = OrnsteinUhlenbeckActionNoise(varargin)
            parser = inputParser;
            addParameter(parser,'InitialAction',0)
            addParameter(parser,'Mean',0)
            addParameter(parser,'MeanAttractionConstant',0.15)
            addParameter(parser,'Variance',0.3)
            addParameter(parser,'VarianceMin',0)  % default is 0 for backward compatibility
            addParameter(parser,'VarianceDecayRate',0)
            addParameter(parser,'SampleTime',-1)
            parse(parser,varargin{:});
            
            obj.InitialAction = parser.Results.InitialAction;
            obj.Mean = parser.Results.Mean;
            obj.MeanAttractionConstant = parser.Results.MeanAttractionConstant;
            obj.Variance = parser.Results.Variance;
            obj.VarianceMin = parser.Results.VarianceMin;
            obj.VarianceDecayRate = parser.Results.VarianceDecayRate;
            obj.SampleTime = parser.Results.SampleTime;
        end
        function obj = set.InitialAction(obj,Value)
            val = @(Value) validateattributes(Value,{'numeric'},{'nonempty','finite'},'','InitialAction');
            checkParam(obj,'InitialAction',val,Value);
            obj.InitialAction = Value;
        end
        function obj = set.Mean(obj,Value)
            val = @(Value) validateattributes(Value,{'numeric'},{'nonempty','finite'},'','Mean');
            checkParam(obj,'Mean',val,Value);
            obj.Mean = Value;
        end
        function obj = set.MeanAttractionConstant(obj,Value)
            val = @(Value) validateattributes(Value,{'numeric'},{'real','positive','nonempty','finite'},'','MeanAttractionConstant');
            checkParam(obj,'MeanAttractionConstant',val,Value);
            obj.MeanAttractionConstant = Value;
        end
        function obj = set.Variance(obj,Value)
            val = @(Value) validateattributes(Value,{'numeric','cell'},{'real','nonnegative','nonempty','finite','nonnan'},'','Variance');
            checkParam(obj,'Variance',val,Value);
            if ~isempty(obj.Variance) || ~isempty(obj.VarianceMin)
                checkParameterValues(obj,Value,'<','VarianceMin','warning','rl:agent:warnNoiseVarianceLtVarianceMin');
            end
            obj.Variance_ = Value;
        end
        function obj = set.VarianceMin(obj,Value)
            val = @(Value) validateattributes(Value,{'numeric','cell'},{'real','nonnegative','nonempty','finite','nonnan'},'','VarianceMin');
            checkParam(obj,'VarianceMin',val,Value);
            if ~isempty(obj.Variance) || ~isempty(obj.VarianceMin)
                checkParameterValues(obj,Value,'>','Variance','warning','rl:agent:warnNoiseVarianceLtVarianceMin');
            end
            obj.VarianceMin_ = Value;
        end
        function obj = set.VarianceDecayRate(obj,Value)
            val = @(Value) validateattributes(Value,{'numeric'},{'real','nonnegative','<=',1,'nonempty'},'','VarianceDecayRate');
            checkParam(obj,'VarianceDecayRate',val,Value);
            obj.VarianceDecayRate = Value;
        end
        function obj = set.SampleTime(obj,Value)
            val = @(Value) validateattributes(Value,{'numeric'},{'real','scalar','nonempty','finite'},'','SampleTime');
            checkParam(obj,'SampleTime',val,Value);
            obj.SampleTime = Value;
        end
        function Options = get.Variance(this)
            Options = this.Variance_;
        end
        function Options = get.VarianceMin(this)
            Options = this.VarianceMin_;
        end
    end
    
    methods (Access = protected)
        function noiseProperties = getNoiseProperties(~)
            noiseProperties = {...
                    'InitialAction',...
                    'Mean',...
                    'MeanAttractionConstant',...
                    'Variance',...
                    'VarianceDecayRate',...
                    'SampleTime'...
                    };
        end
    end
end