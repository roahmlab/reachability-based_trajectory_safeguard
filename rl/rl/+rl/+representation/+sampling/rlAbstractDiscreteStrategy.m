classdef rlAbstractDiscreteStrategy < rl.representation.sampling.rlAbstractSamplingStrategy
    % Discrete uniform sampling object
    
    % Copyright 2019 The MathWorks, Inc.
    
    methods
        function this = rlAbstractDiscreteStrategy(SamplingSpec)
            
            this = this@rl.representation.sampling.rlAbstractSamplingStrategy(SamplingSpec);
            if rl.util.isaSpecType(SamplingSpec, 'continuous')
                error(message('rl:agent:errSamplingStrategySpecMismatch'));
            end
        end
    end
    
    methods (Access = protected)
        function CategoricalProbability = distributionFcn(~, CategoricalProbability)
            % No-op
        end
    end
end

