classdef rlDiscreteSamplingStrategy < rl.representation.sampling.rlAbstractDiscreteStrategy
    % Discrete uniform sampling object
    
    % Copyright 2019 The MathWorks, Inc.
    
    methods
        function this = rlDiscreteSamplingStrategy(SamplingSpec)
            
            this = this@rl.representation.sampling.rlAbstractDiscreteStrategy(SamplingSpec);
        end
    end
    
    methods (Access = protected)
        function Action = sampleImpl(this, Probability, BatchSize, SequenceLength)
            % Return a value sampled from probability vector prob
            % Use in getAction of stochastic discrete policy
            % Probability: 1x1 cell contains 
            % NumProb x BatchSize x SequenceLength
            
            % REVISIT: vectorize to handle batch dim
            if SequenceLength > 1 && BatchSize > 1
                % REVISIT: support batch sequence sampling
                error(message('rl:agent:errDiscreteSamplingBatchSequence'))
            end
            NumExperience = BatchSize * SequenceLength;
            Probability = Probability{1};
            ActionIndex = zeros(1,NumExperience);
            if isa(Probability,'dlarray')
                Probability = extractdata(Probability);
            end
            for ct = 1:NumExperience
                p = Probability(:,ct)' / sum(Probability(:,ct));
                edges = min([0 cumsum(p)],1); % protect against accumulated round-off
                edges(end) = 1;               % get the upper edge exact
                [~, ~, ActionIndex(ct)] = histcounts(rand(),edges);
            end
            % map sampled index to action
            Action = getElementValue(this.SamplingSpec,ActionIndex);
            if ~iscell(Action)
                Action = {Action};
            end
        end
        
        function Model = setModelOutputTypeImpl(~, Model)
            %Model = SETMODELOUTPUTTYPEIMPL(SamplingStrategy, Model) set/
            % validate the model output type to be compatible with the
            % sampling strategy.
            % discrete sampling strategy will add softmax to the end
            % of the model output (output probability)
            
            PostProcessFcn = @(ModelOutput) softmax(ModelOutput);
            Model = setOutputType(Model, 'probability', PostProcessFcn);
        end
    end
end

