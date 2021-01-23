classdef GaussianActionNoise < matlab.mixin.Copyable
    % GAUSSIANACTIONNOISE: Adds Gaussian noise to continuous action
    
    % Copyright 2019 The MathWorks Inc.    
    
    properties
        % Mean expected value of the random variable
        Mean
        
        % Volatility of the random variable
        Variance
        
        % Decay rate of the variance (0 -> no decay, 1-> immediate decay)
        VarianceDecayRate
        
        % Minimum value of variance
        VarianceMin
        
        % Bound of noise samples
        LowerLimit
        UpperLimit
        
        % initial setting of variance for reset
        InitialVariance
    end
    methods        
        function this = GaussianActionNoise(ActionSize, NoiseOpts)
            % build the noise model with the size of the action space, the
            % agent sample time, and the agent noise options
            
            validateattributes(NoiseOpts, {'rl.option.GaussianActionNoise'}, {'scalar', 'nonempty'}, mfilename, 'NoiseOpts', 2);
            
            if ~iscell(ActionSize)
                ActionSize = {ActionSize};
            end
            for i = 1:numel(ActionSize)
                as = ActionSize{i};
                ONES = ones(as);
                try
                    this.Mean{i}              = iCellExtract(NoiseOpts.Mean,i).*ONES;
                    this.Variance{i}          = iCellExtract(NoiseOpts.Variance,i).*ONES;
                    this.VarianceDecayRate{i} = iCellExtract(NoiseOpts.VarianceDecayRate,i).*ONES;
                    this.VarianceMin{i}       = iCellExtract(NoiseOpts.VarianceMin,i).*ONES;
                    this.LowerLimit{i}        = iCellExtract(NoiseOpts.LowerLimit,i).*ONES;
                    this.UpperLimit{i}        = iCellExtract(NoiseOpts.UpperLimit,i).*ONES;
                    this.InitialVariance{i}   = this.Variance{i};
                catch
                    % REVISIT: rename error message?
                    error(message('rl:agent:errOUNoiseInconsistentWithAction',mat2str(ActionSize)));
                end
            end
        end
        
        function action = applyNoise(this,action)
            % apply noise to the action
            
            if iscell(action)
                for i = 1:numel(action)
                    noise = this.Mean{i} + rand(size(action{i})) .* this.Variance{i};
                    noise = min(max(noise,this.LowerLimit{i}),this.UpperLimit{i});
                    action{i} = action{i} + noise;
                end
            else
                noise = this.Mean{1} + rand(size(action)) .* this.Variance{1};
                noise = min(max(noise,this.LowerLimit{1}),this.UpperLimit{1});
                action = action + noise;
            end
            
            % update noise every time step
            update(this)
        end
        
        function reset(this)
            % Resets initial conditions of noise process
            this.Variance = this.InitialVariance;
        end
        
        function update(this)
            % Update noise every time step
            
            for ct = 1:numel(this.Variance)
                if this.VarianceDecayRate{ct}
                    decayedVariance = this.Variance{ct} .* (1-this.VarianceDecayRate{ct});
                    this.Variance{ct} = max(decayedVariance, this.VarianceMin{ct});
                end
            end
        end
    end
end
function val = iCellExtract(val,idx)
    if iscell(val)
        val = val{idx};
    end   
end