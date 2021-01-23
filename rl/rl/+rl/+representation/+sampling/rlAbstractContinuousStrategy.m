classdef rlAbstractContinuousStrategy < rl.representation.sampling.rlAbstractSamplingStrategy
    % Abstract continuous sampling object
    
    % Copyright 2019 The MathWorks, Inc.
    
    methods
        function this = rlAbstractContinuousStrategy(SamplingSpec)
            
            this = this@rl.representation.sampling.rlAbstractSamplingStrategy(SamplingSpec);
            if rl.util.isaSpecType(SamplingSpec, 'discrete')
                error(message('rl:agent:errSamplingStrategySpecMismatch'));
            end
        end
    end
end