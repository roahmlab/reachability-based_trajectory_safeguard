classdef rlPPOAgent < rl.agent.AbstractOnPolicyPGAgent
    % rlPPOAgent: Implements proximal policy optimization agent
    
    % Copyright 2019 The MathWorks Inc.
    
    properties (Access = private)
        % Util to update step tracker
        LastUpdateStep
    end
    
    properties (Constant, Access = private)
        % version indicator for backward compatibility
        Version = 2
    end
    
    methods
        function this = rlPPOAgent(Actor, Critic, Options)
            % Constructor
            
            this = this@rl.agent.AbstractOnPolicyPGAgent(Actor, Critic, Options);
            
            % construct experience buffer and initialize step tracker
            resetImpl(this)
        end
    end
    
    methods(Access = protected)
        function Action = learn(this, Exp)
            % Learn from the current set of experiences where
            % exp = {state,action,reward,nextstate,isdone}
            % return action with exploration.
            
            % store experiences
            appendExperience(this, Exp);
            % update every n steps or if episode terminates
            StepCount = this.EpisodeInfo.StepsTaken;
            BufferLength = StepCount - this.LastUpdateStep + 1;
            DoUpdate = (Exp{5} || BufferLength >= this.AgentOptions.ExperienceHorizon);
            
            if DoUpdate
                % REVISIT: support parallel RNN training
                if hasState(this.Actor)
                    CurrentActorState = getState(this.Actor);
                    this.Actor = resetState(this.Actor);
                end
                
                % compute advantage, critic targets from raw experience
                ProcessedExperience = preprocessExperience(this);
                % train representations over multiple epochs with minibatch 
                % update from processed experience
                learnFromExperiences(this, ProcessedExperience);
                
                if hasState(this)
                    % after the trainig, recover the saved hidden state.
                    this.Actor  = setState(this.Actor, CurrentActorState);
                end
            end
            
            if Exp{5}
                % reset step tracker
                this.LastUpdateStep = 1;
            end
            
            % compute action from the current policy
            % exp is {observation,action,reward,nextobservation,isdone}
            Action = getActionWithExploration(this, Exp{4});
        end
        
        function Actor = setActorLoss(this, Actor)
            if all(arrayfun(@(x) isa(x,'rl.util.rlNumericSpec'), this.ActionInfo))
                Actor = setLoss(Actor,"rl.loss.ppoClipContinuous");
            else
                Actor = setLoss(Actor,"rl.loss.ppoClipDiscrete");
            end
        end
        
        function validateOptionImpl(this,NewOptions)
            validateattributes(NewOptions,{'rl.option.rlPPOAgentOptions'},{'scalar'},'','AgentOptions');
            
            % SequenceLength must be greater than 1 for RNN (PPO treats
            % MiniBatchSize as SequenceLength for RNN)
            if ~isempty(this.Actor) && hasState(this) && NewOptions.MiniBatchSize < 2
                error(message('rl:agent:errPPONotSupportSequenceLengthEq1'))
            end
        end
        
        function trainingOptions = validateAgentTrainingCompatibilityImpl(this,trainingOptions)
            % Validate PPO agent training options compatibility
            
            if ~strcmpi(trainingOptions.Parallelization,'none')
                parallelOptions = trainingOptions.ParallelizationOptions;
                % PPO agent only support send experiences for parallel
                if ~strcmpi(parallelOptions.DataToSendFromWorkers,'Experiences')
                    error(message('rl:general:errParallelSendGradNotSupport'));
                end
                if (parallelOptions.StepsUntilDataIsSent >= trainingOptions.MaxStepsPerEpisode)
                    warning(message('rl:agent:warnPPOStepSendDataFullEpi'));
                    trainingOptions.ParallelizationOptions.StepsUntilDataIsSent = -1;
                end
                if (parallelOptions.StepsUntilDataIsSent ~= this.AgentOptions.ExperienceHorizon)
                    warning(message('rl:agent:warnPPOExpHorizonNotUseInParallel'));
                end
            end
        end
        
        function resetImpl(this)
            % Construct experience buffer and initialize step tracker
            
            this.ExperienceBuffer = rl.util.ExperienceBuffer(this.AgentOptions.ExperienceHorizon, this.ObservationInfo, this.ActionInfo);
            this.ExperienceBuffer.DoValidate = false;
            this.LastUpdateStep = 1;
            
            if ~isempty(this.Actor) && hasState(this)
                % SequenceLength must be greater than 1 for RNN
                % NOTE: repeat check to catch issue when representation is
                % set after construction
                if this.AgentOptions.MiniBatchSize < 2
                    error(message('rl:agent:errPPONotSupportSequenceLengthEq1'))
                end
                this.Actor = resetState(this.Actor);
            end
            if ~isempty(this.Critic) && hasState(this.Critic)
                this.Critic = resetState(this.Critic);
            end
        end
    end
    
    methods (Hidden)
        function ProcessedExperience = preprocessExperience(this, Experience)
            % Compute advantage, loss variables from the experience buffer
            % or from input experience.
            
            if nargin < 2
                % If experience input is not specified, process data
                % already in the experience buffer (series train)
                Experience = [];
            end
            if ~isempty(Experience)
                % If experience input is specified, append it to the 
                % agent's experience buffer (parallel train)
                append(this.ExperienceBuffer, Experience);
            end
            
            % Compute advantage and critic targets
            DiscountFactor = this.AgentOptions.DiscountFactor;
            switch this.AgentOptions.AdvantageEstimateMethod
                case "gae"
                    [Advantages, CriticTargets, BatchExperience] = computeGeneralizedAdvantage(this.ExperienceBuffer, ...
                        this.Critic, DiscountFactor, this.AgentOptions.GAEFactor);
                case "finite-horizon"
                    [Advantages, CriticTargets, BatchExperience] = computeFiniteHorizonAdvantage(this.ExperienceBuffer, ...
                        this.Critic, DiscountFactor);
            end
            % Unpack experience
            Observation = BatchExperience{1};
            Action      = BatchExperience{2};
            % Compute old action probabilities (piOld)
            OldActionProb = evaluate(this.Actor, Observation);
            OldActionProb = OldActionProb{1};
            if isa(this.ActionInfo,'rl.util.rlNumericSpec')
                % continuous action
                % REVISIT: support multi continuous action channels
                % compute piOld(a|mu,sigma)
                OldActionProb = evaluate(this.Actor.SamplingStrategy, OldActionProb, Action{1});
                % loss required old action
                ProcessedExperience.Action = Action{1};
            else
                % discrete action 
                %   pi old is output from actor representation
                %   loss required action index matrix (logical of which action was choose for each obs)
                % REVISIT: support multi discrete action channels
                ActionIndicationMatrix = getElementIndicationMatrix(this.ActionInfo, Action, numel(CriticTargets));
                ProcessedExperience.Action = ActionIndicationMatrix;
            end
            
            % Package processed experience before sending to learner
            ProcessedExperience.Observation = Observation;
            ProcessedExperience.CriticTargets = CriticTargets;
            ProcessedExperience.Advantages = Advantages;
            ProcessedExperience.OldActionProb = OldActionProb;
            
            % Empty experience buffer at the end of the update interval
            reset(this.ExperienceBuffer);
            this.LastUpdateStep = this.EpisodeInfo.StepsTaken + 1;
        end
        
        function learnFromExperiences(this, ProcessedExperience)
            % Train representations over multiple epochs with minibatch 
            % update from processed experience
            
            % Unpack network inputs, loss variables
            Observation = ProcessedExperience.Observation;
            CriticTargets = ProcessedExperience.CriticTargets;
            Advantages = ProcessedExperience.Advantages;
            Action = ProcessedExperience.Action;
            OldActionProb = ProcessedExperience.OldActionProb;
            % NOTE: RNN training only supports single batch
            BufferLength = numel(CriticTargets);
            
            LossVariable.SamplingStrategy = this.Actor.SamplingStrategy;
            LossVariable.EntropyLossWeight = this.AgentOptions.EntropyLossWeight;
            LossVariable.ClipFactor = this.AgentOptions.ClipFactor;
            % get batch dimension or sequence dimension to slice, only
            % applicable if BufferLength > 1
            if hasState(this.Actor)
                ObsDimsToSlice = cellfun(@(x) numel(x) + 2, this.ExperienceBuffer.ObservationDimension','UniformOutput',false);
            else
                ObsDimsToSlice = cellfun(@(x) numel(x) + 1, this.ExperienceBuffer.ObservationDimension','UniformOutput',false);
            end
            for epoch = 1:this.AgentOptions.NumEpoch
                if BufferLength > 1
                    % shuffle dataset if not RNN
                    DoShuffle = ~hasState(this);
                    MiniBatchIdx = rl.internal.dataTransformation.getMiniBatchIdx(BufferLength, this.AgentOptions.MiniBatchSize, DoShuffle);
                    for ct = 1:numel(MiniBatchIdx)
                        % Slice mini batch data
                        SingleBatchIdx = MiniBatchIdx{ct};
                        MiniBatchObs           = rl.internal.dataTransformation.generalSubref(Observation, SingleBatchIdx, ObsDimsToSlice);
                        MiniBatchCriticTargets = rl.internal.dataTransformation.generalSubref(CriticTargets, SingleBatchIdx, ndims(CriticTargets));
                        % REVISIT: support single action channel
                        LossVariable.Action    = rl.internal.dataTransformation.generalSubref(Action, SingleBatchIdx, ndims(Action));
                        LossVariable.OldPolicy = rl.internal.dataTransformation.generalSubref(OldActionProb, SingleBatchIdx, ndims(OldActionProb));
                        LossVariable.Advantage = rl.internal.dataTransformation.generalSubref(Advantages, SingleBatchIdx, ndims(Advantages));
                        
                        % Scale the gradient based on ratio of current minibatch size over specified minibatch size
                        GradScale = single(numel(SingleBatchIdx)/this.AgentOptions.MiniBatchSize);
                        
                        % Update Critic
                        GradVal = gradient(this.Critic, 'loss-parameters', MiniBatchObs, MiniBatchCriticTargets);
                        GradVal = rl.internal.dataTransformation.scaleLearnables(GradVal, GradScale);
                        this.Critic = optimize(this.Critic, GradVal);
                        
                        % Update Actor
                        GradVal = gradient(this.Actor,'loss-parameters',MiniBatchObs,LossVariable);
                        GradVal = rl.internal.dataTransformation.scaleLearnables(GradVal, GradScale);
                        this.Actor = optimize(this.Actor, GradVal);
                    end
                else
                    % Avoid minibatch slicing if only 1 experience
                    LossVariable.Action    = Action;
                    LossVariable.OldPolicy = OldActionProb;
                    LossVariable.Advantage = Advantages;
                        
                    % Scale the gradient based on ratio of current minibatch size over specified minibatch size
                    GradScale = single(1/this.AgentOptions.MiniBatchSize);
                    
                    % Update Critic
                    GradVal = gradient(this.Critic, 'loss-parameters', Observation, CriticTargets);
                    GradVal = rl.internal.dataTransformation.scaleLearnables(GradVal, GradScale);
                    this.Critic = optimize(this.Critic, GradVal);
                    
                    % Update Actor
                    GradVal = gradient(this.Actor,'loss-parameters',Observation, LossVariable);
                    GradVal = rl.internal.dataTransformation.scaleLearnables(GradVal, GradScale);
                    this.Actor = optimize(this.Actor, GradVal);
                end
            end
        end
    end
    
    methods (Static)
        function obj = loadobj(s)
            if isstruct(s)
                if ~isfield(s,'Version')
                    % version 1 does not have Version field
                    % In version 2,
                    %   - Critic changes from rlRepresentation to rlValueRepresentation
                    %   - Actor changes from rlRepresentation to rlStochasticActorRepresentation
                    %   - ExperienceBuffer requires obs and act dims inputs
                    %   but will always reconstruct since agent is on-policy
                    [Actor, Critic] = rl.agent.AbstractOnPolicyPGAgent.redirectV1Rep(s);
                    obj = rl.agent.rlPPOAgent(Actor,Critic,s.AgentOptions_);
                end
            else
                obj = s;
            end
        end
    end
end