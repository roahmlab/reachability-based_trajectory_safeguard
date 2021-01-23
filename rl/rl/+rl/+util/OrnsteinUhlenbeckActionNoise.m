classdef OrnsteinUhlenbeckActionNoise < matlab.mixin.Copyable
    % ORNSTEINUHLENBECKACTIONNOISE: Adds noise to continuous random
    % variables with mean reverting properties
    
    % Copyright 2017-2018 The MathWorks Inc.    
    
    properties        
        % How quickly random variable reverts to the mean
        MeanAttractionConstant
        
        % Mean expected value of the random variable
        Mean
        
        % Volatility of the random variable
        Variance
        
        % Minimum value of variance
        VarianceMin
        
        % decay rate of the variance (0 -> no decay, 1-> immediate decay)
        VarianceDecayRate
        
        % Time elapsed between samples
        SampleTime
        
        % Previous state of the process
        Xprev
        
        % Initial state of the process
        X0
        
        % initial setting of variance for reset
        InitialVariance
    end
    methods        
        function this = OrnsteinUhlenbeckActionNoise(actionSize,noiseOpts,agentTs)
            % build the noise model with the size of the action space, the
            % agent sample time, and the agent noise options
            
            ts = noiseOpts.SampleTime;
            if ts == -1
                ts = agentTs;
            end
            
            if ~iscell(actionSize)
                actionSize = {actionSize};
            end
            for i = 1:numel(actionSize)
                as = actionSize{i};
                ONES = ones(as);
                try
                    this.MeanAttractionConstant{i} = localCellExtract(noiseOpts.MeanAttractionConstant,i).*ONES;
                    this.Mean{i}                   = localCellExtract(noiseOpts.Mean,i).*ONES;
                    this.Variance{i}               = localCellExtract(noiseOpts.Variance,i).*ONES;
                    this.VarianceMin{i}            = localCellExtract(noiseOpts.VarianceMin,i).*ONES;
                    this.VarianceDecayRate{i}      = localCellExtract(noiseOpts.VarianceDecayRate,i).*ONES;
                    this.SampleTime{i}             = localCellExtract(ts,i).*ONES;
                    this.X0{i}                     = localCellExtract(noiseOpts.InitialAction,i).*ONES;
                    this.InitialVariance{i}        = this.Variance{i};
                catch
                    error(message('rl:agent:errOUNoiseInconsistentWithAction',mat2str(actionSize)));
                end
            end
            this.reset();
        end
        
        function action = applyNoise(this,action)
            % apply noise to the action
            noise = update(this);
            if iscell(action)
                for i = 1:numel(action)
                    action{i} = action{i} + noise{i};
                end
            else
                action = action + noise{1};
            end
        end
        
        function x = update(this)
            % Update noise every time step
            for i = numel(this.Mean):-1:1
                x{i} = this.Xprev{i} + ...
                    this.MeanAttractionConstant{i}.*(this.Mean{i} - this.Xprev{i}).*this.SampleTime{i} + ...
                    this.Variance{i}.*randn(size(this.Mean{i})).*sqrt(this.SampleTime{i});
                % update the last state
                this.Xprev{i} = x{i};
                % update the variance
                decayedVariance = this.Variance{i} .* (1-this.VarianceDecayRate{i});
                this.Variance{i} = max(decayedVariance, this.VarianceMin{i});
            end
        end
        
        function reset(this)
            % Resets initial conditions of noise process
            this.Xprev = this.X0;
            this.Variance = this.InitialVariance;
        end
    end
    
    methods(Static)
       function obj = loadobj(s)
          if isstruct(s)
              obj = rl.util.OrnsteinUhlenbeckActionNoise;
              obj.MeanAttractionConstant = s.MeanAttractionConstant;
              obj.Mean = s.Mean;
              obj.Variance = s.Variance;
              obj.VarianceDecayRate = s.VarianceDecayRate;
              obj.SampleTime = s.SampleTime;
              obj.Xprev = s.Xprev;
              obj.X0 = s.X0;
              obj.InitialVariance = s.InitialVariance;
          else
              obj = s;
          end
          if isempty(s.VarianceMin)
              % set default to 0 if VarianceMin is not present
              % VarianceMin has same dimension with Variance
              obj.VarianceMin = cellfun(@(x) x*0,obj.Variance,'UniformOutput',false);
          else
              obj.VarianceMin = s.VarianceMin;
          end
       end
    end
end
function val = localCellExtract(val,idx)
    if iscell(val)
        val = val{idx};
    end   
end