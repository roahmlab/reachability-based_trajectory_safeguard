classdef ExperienceBuffer < matlab.mixin.Copyable
    % EXPERIENCEBUFFER: Stores agent experiences (states, actions,
    % and rewards) in a buffer. Agents sample data from this buffer
    % during training.
    %
    % Examples:
    %
    % 1) Append experiences to ExperienceBuffer.
    %    buffer = rl.util.ExperienceBuffer(5,{[1 1]},{[1 1]}); % create buffer of capacity 5
    %    exp1 = {{1},{1},1,{2},0}; % experience to be appended
    %    append(buffer,{exp1}); % append experience to buffer
    %    exps = {{{2},{2},-1,{3},0},{{3},{3},2,{4},0},{{4},{4},-2,{5},1}};
    %    append(buffer,exps); % append experiences to buffer
    %
    % 2) Create mini-batch of sampled experiences with n-step return.
    %    miniBatch = createSampledExperienceMiniBatch(ExperienceBuffer,...
    %                                                 BatchSize,...
    %                                                 DiscountFactor,...
    %                                                 NStepsToLookAhead);
    %    miniBatch = createSampledExperienceMiniBatch(buffer,3);
    %    miniBatch = createSampledExperienceMiniBatch(buffer,3,0.8,2);
    % 
    % 3) Create mini-batch of experiences with returns.
    %    miniBatch = createExperienceWithReturnMiniBatch(ExperienceBuffer,...
    %                                                    DiscountFactor,...
    %                                                    NumOfEpisodes);
    %    miniBatch = getExperienceWithReturn(buffer,0.9);
    %    miniBatch = getExperienceWithReturn(buffer,0.9,1);
    %
    % 4) Create mini-batch of indexed experiences with n-step return.
    %    miniBatch = createNStepExperienceMiniBatch(ExperienceBuffer,...
    %                                                 Indexes,...
    %                                                 DiscountFactor,...
    %                                                 NStepsToLookAhead);
    %    miniBatch = createNStepExperienceMiniBatch(buffer,[1 2],0.8,2);
    %
    % 5) Reset ExperienceBuffer.
    %    reset(buffer);
    %
    % 6) Compute sampled experiences.
    %    expNew = getSampledExperience(ExperienceBuffer,...
    %                                 BatchSize,...
    %                                 DiscountFactor,...
    %                                 NStepsToLookAhead);
    %    expNew = getSampledExperience(buffer,3);
    %    expNew = getSampledExperience(buffer,3,0.8,2);
    %
    % 7) Compute experiences with returns.
    %    expNew = getExperienceWithReturn(ExperienceBuffer,...
    %                                     DiscountFactor,...
    %                                     NumOfEpisodes);
    %    expNew = getExperienceWithReturn(buffer,0.9);
    %    expNew = getExperienceWithReturn(buffer,0.9,1);
    % 
    %  8) Stack experiences into mini-batch.
    %     miniBatch = buffer.batchExperience(expNew);
    % 
    % ExperienceBuffer objects are created and managed automatically for
    % standard Reinforcement Learning Toolbox agents. You can also create
    % and manage the experience buffer in your own custom agents.
    
    % Copyright 2017-2018 The MathWorks Inc.
    
    properties
        % Validate input arguments to be true or false
        DoValidate = true
    end
    
    properties (SetAccess = private)
        % Experience buffer Size
        Capacity {mustBeNumeric, mustBePositive, mustBeReal}
        
        % Index to next location to fill
        NextIndex = 1
        
        % Observation and action dimension (for validate and sample data)
        ObservationDimension
        ActionDimension
    end
    
    properties(Dependent)
        % Current number of elements in the circular buffer
        Length
    end
    
    properties(Access = protected)
        % Actual storage for the circular buffer
        Memory = {}
        % Flag for is full
        IsFull = false;
        % save memory property with the buffer
        SaveMemoryWithBuffer (1,1) logical = false
    end
    
    properties(Access = private)
        % Sparse isdone indicator of size 1 x obj.Length to improve data
        % sampling performance (mostly sequence data)
        % NOTE: get update in append_
        IsdoneVector
    end
    
    methods
        %%
        function obj = ExperienceBuffer(Capacity, ObservationDimension, ActionDimension)
            % Constructor
            
            validateattributes(Capacity,{'numeric'},{'scalar','integer','positive','finite'},'','Capacity');
            obj.Capacity = Capacity;
            obj.Memory = cell(1, Capacity);
            obj.NextIndex = 1;
            obj.IsFull = false;
            obj.IsdoneVector = sparse(1, Capacity);
            
            % Take empty ObservationDimension and ActionDimension to avoid 
            % error when loading empty buffer from V1 without observation 
            % and action dimension (backward compatibility)
            if ~isempty(ObservationDimension) && ~isempty(ActionDimension)
                validateattributes(ObservationDimension,{'cell','rl.util.RLDataSpec'},{'vector','nonempty'},'','ObservationDimension');
                validateattributes(ActionDimension,{'cell','rl.util.RLDataSpec'},{'vector','nonempty'},'','ActionDimension');
                
                % parse data dimension spec
                if isa(ObservationDimension, 'rl.util.RLDataSpec')
                    ObservationDimension = {ObservationDimension.Dimension};
                end
                if isa(ActionDimension, 'rl.util.RLDataSpec')
                    ActionDimension = {ActionDimension.Dimension};
                end
                ObservationDimension = iCellify(ObservationDimension);
                ActionDimension = iCellify(ActionDimension);
                
                % validate each entry of dimension cell are non scalar
                validateFcn = @(x) validateattributes(x, {'numeric'}, {'nonsparse', 'row', 'finite'}, mfilename, 'ObservationDimension', 2);
                cellfun(validateFcn, ObservationDimension)
                validateFcn = @(x) validateattributes(x, {'numeric'}, {'nonsparse', 'row', 'finite'}, mfilename, 'ActionDimension', 3);
                cellfun(validateFcn, ActionDimension)
                
                % force column cell array
                obj.ObservationDimension = ObservationDimension(:);
                obj.ActionDimension      = ActionDimension(:);
            else
                obj.DoValidate = true;
            end
        end
        %%
        function reset(obj)
            obj.Memory = cell(1, obj.Capacity);
            obj.NextIndex = 1;
            obj.IsFull = false;
            obj.IsdoneVector = sparse(1, obj.Capacity);
        end
        %%
        function s = saveobj(obj)
            % modify the save process depending if SaveMemoryWithBuffer is
            % off
            s.Capacity = obj.Capacity;
            s.DoValidate = obj.DoValidate;
            s.SaveMemoryWithBuffer = obj.SaveMemoryWithBuffer;
            s.ObservationDimension = obj.ObservationDimension;
            s.ActionDimension = obj.ActionDimension;
            if obj.SaveMemoryWithBuffer
                s.NextIndex = obj.NextIndex;
                s.IsFull = obj.IsFull;
                s.Memory = obj.Memory;
                s.IsdoneVector = obj.IsdoneVector;
            else
                s.NextIndex = 1;
                s.IsFull = false;
                s.Memory = {};
                s.IsdoneVector = sparse(1, s.Capacity);
            end
            s.IsFull = obj.IsFull;
            % version indicator
            s.Version = 2;
        end
        %%
        function append(obj,DataArray)
            % REVISIT: Support append batch data (not array of single data)
            
            if isempty(obj.ObservationDimension) && isempty(obj.ActionDimension) && obj.Length == 1
                % if buffer is constructed without obs and act dim, extract
                % the information from the first experience to the buffer
                validateattributes(DataArray,{'cell'},{'nonempty'},'','Experience');
                [obj.ObservationDimension, obj.ActionDimension] = iGetDimInfoFromExp(DataArray{1});
                DataArray{1} = iCellifyExperience(DataArray{1});
            end
            
            if obj.DoValidate
                validateattributes(DataArray,{'cell'},{'nonempty'},'','Experience');
                for ct = 1:length(DataArray)
                    validateExperience(obj,DataArray{ct});
                end
            end
            append_(obj, DataArray);
        end
        %%
        function minibatch = createSampledExperienceMiniBatch(obj,BatchSize,varargin)
            % create minibatch of sampled experiences
            dataArray = getSampledExperience(obj,BatchSize,varargin{:});
            minibatch = batchExperience(obj,dataArray);
        end
        %%
        function [BatchSequenceData,BatchSequenceMask] = createSampledExperienceMiniBatchSequence(obj,BatchSize,SequenceLength)
            % return BatchSize*SequenceLength of experience and mask (see below). From the
            % buffer, sample randomly BatchSize episodes. From each
            % episode, sample randomly SequenceLength continuous exp
            %
            % BatchSequenceMask: [1 × Batch size × sequenceLength] logic
            % Each sampled trajectory in a batch has a different length. We pad inputs
            % if the trajectory is shorter than specified the trajectory
            % length. If an element of BatchSequenceMask is true, the corresponding
            % input is the original input (non-padded input). If it is false,
            % the input is a padded input. The padded inputs
            % are ignored in the gradient computation.
            
            % REVISIT: can improve performance if we return cell of 
            % 1xBatchSizexSequenceLength instead of 1x(BatchSize*SequenceLength)
            if obj.DoValidate
                narginchk(3,3);
                validateattributes(BatchSize,{'numeric'},{'scalar','integer','positive'},'','BatchSize');
                validateattributes(SequenceLength,{'numeric'},{'scalar','integer','positive','>',1},'','SequenceLength');                
            end
            
            % get all terminal and starting idexes from the buffer
            if obj.IsFull
                % change the index system to represent a FIFO buffer
                % NewIdxSystem = [OldestDataIndex:obj.Capacity (1:NewestDataIndex)+obj.Capacity];
                [NewestDataIndex,OldestDataIndex] = findNewOldDataIdxes(obj);
                ShiftedIsdoneVector = [obj.IsdoneVector(OldestDataIndex:end) obj.IsdoneVector(1:NewestDataIndex)];
                TerminalIdx = find(ShiftedIsdoneVector) + OldestDataIndex-1;
                % get all start idexes from the buffer (not include the last unfinished episode)
                StartIdx = [OldestDataIndex TerminalIdx(1:end-1)+1];
            else
                TerminalIdx = findEpisodeIsdoneIdx(obj);
                % get all start idexes from the buffer
                StartIdx = [1 TerminalIdx(1:end-1)+1];
            end
            EpisodeLength = TerminalIdx - StartIdx + 1;
            NumEpisode = length(TerminalIdx);

            % We want at least N trajectories (the number of episodes) where
            % N == Batchsize before it starts training. If the memory size
            % is small, it never gets N. If memory is already full, it
            % starts training even N < Batchsize.
            if NumEpisode < BatchSize && ~obj.IsFull
                BatchSequenceData = {};
                Masks = {};
                NumEpisodesInBatch=BatchSize;
            else
                % sample episodes                
                NumEpisodesInBatch = min(NumEpisode, BatchSize); 
                SampledEpisode = randperm(NumEpisode,NumEpisode) ;                
                Idxes = zeros(1,NumEpisodesInBatch*SequenceLength);
                Masks = true(1, NumEpisodesInBatch*SequenceLength);

                for ct = 1:NumEpisodesInBatch
                    EpisodeNumber = SampledEpisode(ct);
                    % create padding and masking
                    currentEpisodeLength = EpisodeLength(EpisodeNumber);
                    if currentEpisodeLength < SequenceLength
                        % If the sample trajectory (current) in the batch is shorter than the specified length 
                        % (Sequence Length), we need to padd experience. The padded experience
                        % will be ignored in the gradient computation. 

                        % currentStartIdxInIdxes is the index that specified
                        % the beginning of the sample trajectory in
                        % "Idxes" while StartIdx is the index in
                        % the replay buffer.
                        currentStartIdxInIdxes = (ct-1)*SequenceLength+1;

                        % currentEndIdxInIdxes is the index that specified
                        % the end of the sample trajectory in
                        % Idxes.
                        currentEndIdxInIdxes = currentStartIdxInIdxes+currentEpisodeLength-1;

                        if obj.IsFull
                            % mod() will shift idx to original buffer idx system
                            Idxes(currentStartIdxInIdxes:currentEndIdxInIdxes) = mod(StartIdx(EpisodeNumber):StartIdx(EpisodeNumber)+currentEpisodeLength-1,obj.Capacity);
                        else
                            Idxes(currentStartIdxInIdxes:currentEndIdxInIdxes) = StartIdx(EpisodeNumber):StartIdx(EpisodeNumber)+currentEpisodeLength-1;
                        end
                        % use the last Idx for padding (use the input value at the end of trajectory as padding) and use false for
                        % for masking. When mask is true, it indicates the original input (non-padded input).
                        % When mask is false, it indicates a PADDED input. 
                        Idxes(currentEndIdxInIdxes+1:ct*SequenceLength) = Idxes(currentEndIdxInIdxes);
                        Masks(currentEndIdxInIdxes+1:ct*SequenceLength) = false;
                    else
                        % sample a starting index from the choosen episode
                        % (remove last SequenceLength indexes from current
                        % episode indexes)
                        ExperienceIdxSet = (StartIdx(EpisodeNumber):(TerminalIdx(EpisodeNumber)-SequenceLength+1));
                        RandomStartIdx = StartIdx(EpisodeNumber)-1 + randi(numel(ExperienceIdxSet));
                        if obj.IsFull
                            % mod() will shift idx to original buffer idx system
                            Idxes((ct-1)*SequenceLength+1:ct*SequenceLength) = mod(RandomStartIdx:(RandomStartIdx+SequenceLength-1),obj.Capacity);
                        else
                            Idxes((ct-1)*SequenceLength+1:ct*SequenceLength) = RandomStartIdx:(RandomStartIdx+SequenceLength-1);

                        end
                    end
                end

                if obj.IsFull
                    % When you use mod, the last index of the capcity will
                    % be zero. Then, we need to manually change the index 0
                    % to the last index.
                    Idxes(Idxes == 0) = obj.Capacity;
                end            
                BatchSequenceData = obj.Memory(Idxes);

            end
            [BatchSequenceData, BatchSequenceMask] = batchSequenceExperience(obj,BatchSequenceData,NumEpisodesInBatch,SequenceLength,Masks);
        end        
        %%
        function minibatch = createExperienceWithReturnMiniBatch(obj,gamma,varargin)
            % create minibatch of experiences with return
            dataArray = getExperienceWithReturn(obj,gamma,varargin{:});
            minibatch = batchExperience(obj,dataArray);
        end
        %%
        function minibatch = createNStepExperienceMiniBatch(obj,idxes,gamma,n)
            % create minibatch of experiences locatated at idxes in buffer 
            % with n-step look-ahead return
            
            if obj.DoValidate
                validateattributes(gamma,{'numeric'},{'scalar','positive','<=',1,'>',0},'','DiscountFactor');
                validateattributes(n,{'numeric'},{'scalar','integer','positive','finite'},'','N-step look-ahead');
                maxIndex = max(1,obj.Capacity*obj.IsFull + obj.Length*(~obj.IsFull));
                validateattributes(idxes,{'numeric'},{'vector','integer','positive','<=',maxIndex},'','Experience index');
            end
            if obj.Length < 1
                minibatch = [];
            else
                dataArray = convertToNStepExperiences(obj,idxes,gamma,n);
                minibatch = batchExperience(obj,dataArray);
            end
        end
        
        %% Get methods
        function bufferSize = get.Length(obj)
            % Returns the current number of elements in buffer
            if obj.IsFull
                bufferSize = obj.Capacity;
            else
                bufferSize = obj.NextIndex - 1;
            end
        end
        %% Set methods
        function set.DoValidate(obj,Value)
            validateattributes(Value,{'logical','numeric'},{'scalar','real'},'','DoValidate');
            obj.DoValidate = logical(Value);
        end
        
    end
    
    methods (Hidden)
        %%
        function setSaveMemoryWithBuffer(obj,val)
            obj.SaveMemoryWithBuffer = val;
        end
        function val = getSaveMemoryWithBuffer(obj)
            val = obj.SaveMemoryWithBuffer;
        end
        %%
        function data = getMemory(obj)
            % Returns data in Memory
            data = obj.Memory;
        end
        %%
        function data = getLastNData(obj,N)
            % Returns last N data in Memory
            
            N = min(obj.Capacity,N);
            [NewestDataIndex,~] = findNewOldDataIdxes(obj);
            startIdx = NewestDataIndex - N + 1;
            if ~obj.IsFull
                if N >= NewestDataIndex
                    data = obj.Memory(1:NewestDataIndex);
                else
                    data = obj.Memory(startIdx:NewestDataIndex);
                end
            else
                if N <= NewestDataIndex
                    data = obj.Memory(startIdx:NewestDataIndex);
                else
                    data = [obj.Memory(obj.Capacity-abs(startIdx):obj.Capacity) obj.Memory(1:NewestDataIndex)];
                end
            end
        end
        %%
        function resize(obj,newLength)
            % Change length of buffer while keeping old memory
            
            if newLength ~= obj.Capacity
                if newLength > obj.Capacity
                    oldData = getLastNData(obj,obj.Capacity);
                else
                    warning(message('rl:general:warnDataLostResize'));
                    oldData = getLastNData(obj,newLength);
                end
                obj.Capacity =  newLength;
                obj.Memory = cell(1, newLength);
                obj.NextIndex = 1;
                obj.IsFull = false;
                obj.append_(oldData);
            end
        end
        %%
        function data = getIsFull(obj)
            % Returns data in Memory
            data = obj.IsFull;
        end
        %%
        function [dataArray,idxes] = getSampledExperience(obj,BatchSize,varargin)
            % Get sampled experiences
            % dataArray = getSampledExperience(obj,BatchSize,gamma,n);
            if obj.DoValidate
                narginchk(2,4);
                validateattributes(BatchSize,{'numeric'},{'scalar','integer','positive'},'','MiniBatchSize');
                if nargin>2
                    if nargin~=4
                        error((message('rl:general:errInputsNStep')));
                    end
                    gamma = varargin{1};
                    validateattributes(gamma,{'numeric'},{'scalar','positive','<=',1,'>',0},'','DiscountFactor');
                    n = varargin{2};
                    validateattributes(n,{'numeric'},{'scalar','integer','positive','finite'},'','N-step look-ahead');
                end
            end
            [dataArray,idxes] = getSampledExperience_(obj,BatchSize,varargin{:});
        end
        %%
        function dataArray = getExperienceWithReturn(obj,gamma,varargin)
            % Get experiences with return
            % dataArray = getExperienceWithReturn(obj,gamma,n);
            if obj.DoValidate
                narginchk(2,3);
                validateattributes(gamma,{'numeric'},{'scalar','positive','<=',1,'>',0},'','DiscountFactor');
                if nargin>2
                    validateattributes(varargin{1},{'numeric'},{'scalar','integer','positive','finite'},'','Number of episodes');
                end
            end
            dataArray = getExperienceWithReturn_(obj,gamma,varargin{:});
        end
        %%
        function convertedExperiences = convertToNStepExperiences(obj,idxes,gamma,n)
            % Converts sampled raw data in ReplayMemory to computed data
            % with n-step look-ahead.
            NewestDataIndex = findNewOldDataIdxes(obj);% index for newest data in memory
            numSamples = length(idxes);                % number of sampled experiences
            convertedExperiences = cell(1,numSamples); % converted experiences
            % Convert sampled experiences
            for j = 1:numSamples
                % Initialization
                idx = idxes(j);                         % sampled experience index
                rewards = zeros(n,1);                   % vector of rewards for preview buffer
                maxIdx = min(idx+n-1,...                % maximum index for preview buffer
                    NewestDataIndex + (idx>NewestDataIndex)*obj.Capacity);
                % Find index for last state and preview length
                for ct = idx:maxIdx
                    memoryIndex = mod(ct -1, obj.Capacity) + 1;
                    rewards(ct-idx+1) = obj.Memory{memoryIndex}{3};
                    if obj.Memory{memoryIndex}{5}>=1 || (ct == maxIdx)
                        previewLength = ct-idx+1;
                        lastIndex = memoryIndex;
                        break;
                    end
                end
                % Get rewards, weights vector to compute return
                previewRewards = rewards(1:previewLength);
                weights = gamma.^(0:previewLength-1);
                computedReturn = weights*previewRewards;
                % Converted experience with grouped elements
                convertedExperiences{j} = { ...
                    obj.Memory{idx}{1},...
                    obj.Memory{idx}{2},...
                    computedReturn,...
                    obj.Memory{lastIndex}{4},...
                    obj.Memory{lastIndex}{5}};
            end
        end
        %%
        function [Advantage, TDTarget, BatchExperience] = computeFiniteHorizonAdvantage(obj, StateValueEstimator, DiscountFactor)
            % Vectorized finite horizon advantage estimator (A3C)
            % REVISIT: current implementation supports single episode
            
            BatchExperience = getBatchExperience(obj,hasState(StateValueEstimator));
            
            % Unpack experience
            Observation       = BatchExperience{1};
            Reward            = BatchExperience{3};
            NextObservation   = BatchExperience{4};
            IsDone            = BatchExperience{5};
            SequenceLength = numel(Reward);
            
            % Estimate state value
            CurrentStateValue = getValue(StateValueEstimator,Observation);
            if hasState(StateValueEstimator)
                % slice on the sequence dimension
                DimToSlice = cellfun(@(x) numel(x) + 2, obj.ObservationDimension', 'UniformOutput', false);
            else
                % slice on the batch dimension
                DimToSlice = cellfun(@(x) numel(x) + 1, obj.ObservationDimension', 'UniformOutput', false);
            end
            LastObservation = rl.internal.dataTransformation.generalSubref(NextObservation, SequenceLength, DimToSlice);
            
            if IsDone(end) == 1
                TargetStateValue = 0;
            else
                TargetStateValue = getValue(StateValueEstimator, LastObservation);
            end
            
            % Vectorize Finite Horizon
            % Adv    = ReturnNstep + gamma^n*V[t+n] - V[t]
            
            % ReturnNstep(1) = r(1) + gamma*r(2) + gamma^2*r(3) + gamma^3*r(4)
            % ReturnNstep(2) =              r(2) + gamma^1*r(3) + gamma^2*r(4)
            % ReturnNstep(3) =                             r(3) + gamma^1*r(4)
            % ReturnNstep(4) =                                            r(4)
            % ...
            % ReturnNstep    = [r(1) r(2) ... r(4)] * [      1         0         0    0
            %                                          gamma^1         1         0    0
            %                                          gamma^2   gamma^1         1    0
            %                                          gamma^3   gamma^2   gamma^1    1]
            % ReturnNstep    =        Reward        *  DiscountWeights
            
            WeightsMatrix = repmat((0:SequenceLength-1)',1,SequenceLength) - (0:SequenceLength-1);
            %WeightsMatrix =
            %   [0   -1   -2   -3   -4
            %    1    0   -1   -2   -3
            %    2    1    0   -1   -2
            %    3    2    1    0   -1
            %    4    3    2    1    0]
            
            DiscountWeights = tril(DiscountFactor .^ WeightsMatrix);
            % With A = gamma, DiscountWeights =
            %  [  1     0     0    0
            %   A^1     1     0    0
            %   A^2   A^1     1    0
            %   A^3   A^2   A^1    1]
            
            ReturnNstep = Reward(:)' * DiscountWeights;
            
            % Vectorize gamma^n*V[t+n]
            DiscountedTargetStateValue = DiscountFactor.^(SequenceLength:-1:1) .* TargetStateValue;
            % Temporal different target = n-step return + gamma^n*V[t+n]
            TDTarget = ReturnNstep + DiscountedTargetStateValue;
            TDTarget = reshape(TDTarget,size(CurrentStateValue));
            % Advantages = TDTarget - CurrentStateValue
            Advantage = TDTarget - CurrentStateValue;
        end
        %%
        function [Advantage, TDTarget, BatchExperience] = computeGeneralizedAdvantage(obj, StateValueEstimator, DiscountFactor, GAEFactor)
            % Vectorized generalized advantage estimator (GAE)
            % REVISIT: current implementation supports single episode
            
            BatchExperience = getBatchExperience(obj,hasState(StateValueEstimator));
            
            % Unpack experience
            Observation       = BatchExperience{1};
            Reward            = BatchExperience{3};
            NextObservation   = BatchExperience{4};
            IsDone            = BatchExperience{5};
            SequenceLength = numel(Reward);
            
            % Estimate current and next state values
            CurrentStateValue = getValue(StateValueEstimator, Observation);
            NextStateValue = getValue(StateValueEstimator, NextObservation);
            NextStateValue(IsDone == 1) = 0; % early termination
            
            % Vectorized GAE Advantages
            % TDError = [TDError(1) TDError(2) ... TDError(4)]
            TDError = Reward + ...
                reshape(DiscountFactor * NextStateValue - CurrentStateValue, size(Reward));
            if GAEFactor == 0
                % If GAEFactor == 0, similar to 1 step look ahead (or TD0)
                Advantage = TDError;
            else
                % Adv(1) = TDError(1) + A*TDError(2) + A^2*TDError(3) + A^3*TDError(4)
                % Adv(2) =                TDError(2) + A^1*TDError(3) + A^2*TDError(4)
                % Adv(3) =                                 TDError(3) + A^1*TDError(4)
                % Adv(4) =                                                  TDError(4)
                % ...
                % Adv    = [TDError(1) TDError(2) ... TDError(4)] * [  1     0     0    0
                %                                                    A^1     1     0    0
                %                                                    A^2   A^1     1    0
                %                                                    A^3   A^2   A^1    1]
                % Adv    =                   TDError              *  DiscountWeights
                
                WeightsMatrix = repmat((0:SequenceLength-1)',1,SequenceLength) - (0:SequenceLength-1);
                %WeightsMatrix =
                %   [0   -1   -2   -3   -4
                %    1    0   -1   -2   -3
                %    2    1    0   -1   -2
                %    3    2    1    0   -1
                %    4    3    2    1    0]
                DiscountWeights = tril((DiscountFactor*GAEFactor) .^ WeightsMatrix);
                % With A = DiscountFactor*GAELambda, DiscountWeights =
                %  [  1     0     0    0
                %   A^1     1     0    0
                %   A^2   A^1     1    0
                %   A^3   A^2   A^1    1]
                Advantage = TDError(:)' * DiscountWeights;
            end
            % Temporal different target = Advantage[s] + V[s]
            Advantage = reshape(Advantage, size(CurrentStateValue));
            TDTarget = Advantage + CurrentStateValue;
        end
        
        %%
        function [IsdoneIdx,NumEpisodeInBuffer] = findEpisodeIsdoneIdx(obj)
            % Return terminal indexes of each episode in the buffer
            IsdoneIdx = find(obj.IsdoneVector);
            NumEpisodeInBuffer = numel(IsdoneIdx);
        end
        %%
        function [NewestDataIndex,OldestDataIndex] = findNewOldDataIdxes(obj)
            % Find the index for newest data and oldest data in memory
            if obj.IsFull
                OldestDataIndex = obj.NextIndex;
                NewestDataIndex = (obj.NextIndex-1) + (obj.NextIndex == 1) * obj.Capacity;
            else
                OldestDataIndex = 1;
                NewestDataIndex = obj.NextIndex-1;
            end
        end
        
        %%
        function BatchExperience = getBatchExperience(obj, IsSequence)
            % Return 1x5 cell array of obs,action,reward,nextObs,isDone
            % Each cell is a numeric matrix, converted from 1xN cell arrays
            % in Memory 
            % IsSequence indicate the output format
            %   - true:  format [DxB]
            %   - false: format [DxBxT]. B = 1 in this case
            % D: data size
            % B: batch size
            % T: sequence lenght
            
            if nargin < 2
                IsSequence = false;
            end
            
            idx = 1:obj.Length;
            if obj.Length<1
                BatchExperience = [];
            else
                if IsSequence
                    BatchExperience = batchSequenceExperience(obj,obj.Memory(idx),1,obj.Length,{});
                else
                    BatchExperience = batchExperience(obj,obj.Memory(idx));
                end
            end
        end
    end
    
    methods (Access = private)
        %%
        function convertedExperiences = convertToMonteCarloExperiences(obj,idxes,gamma,n)
            % Converts raw data for the latest n episodes in ReplayMemory
            % to computed data with Monte-Carlo returns
            bufferLength = idxes(1)-idxes(n+1);          % length of experiences for n episodes
            convertedExperiences = cell(1,bufferLength); % experiences for n episodes
            offset = idxes(n+1);
            % Every episode (in reverse order)
            for ct = n+1:-1:2
                episodeStart = idxes(ct)+1;     % Start index of the episode
                episodeEnd = idxes(ct-1);       % End index of the episode
                % Index for last state of the episode
                lastIndex = mod(episodeEnd -1, obj.Capacity) + 1;
                % Compute return for the episode
                runningSum = 0;                 % for Monte-Carlo computation
                for idx = episodeEnd:-1:episodeStart
                    % Elements of converted experience
                    memoryIdx = mod(idx -1, obj.Capacity) + 1;
                    runningSum = runningSum * gamma + obj.Memory{memoryIdx}{3};
                    computedReturn = runningSum;
                    % Converted experience with grouped elements
                    convertedExperiences{idx-offset} = {...
                        obj.Memory{memoryIdx}{1},...
                        obj.Memory{memoryIdx}{2},...
                        computedReturn,...
                        obj.Memory{lastIndex}{4},...
                        1};
                end
            end
        end
        %%
        function [idxes,numEpisodeEdges] = findEpisodeIdxes(obj,maxIndex,OldestDataIndex,n)
            % Update Indexes for end of episodes to be learned
            % numEpisode: number of episodes found to be learned
            numEpisodeEdges = 0;
            % idxes: end indexes for each episode to be learned
            idxes = zeros(1,n+1);
            % Find the index when an episode ends
            for ct = maxIndex:-1:OldestDataIndex
                memoryIdx = mod(ct -1, obj.Capacity) + 1;
                if obj.Memory{memoryIdx}{5}>0 || (ct == 1)
                    numEpisodeEdges = numEpisodeEdges + 1;
                    idxes(numEpisodeEdges) = ct - (ct==1);
                    % Add one edge when the first experience is an episode
                    if obj.Memory{memoryIdx}{5}>0 && (ct == 1)
                        idxes(numEpisodeEdges) = 1;
                        numEpisodeEdges = numEpisodeEdges + 1;
                    end
                    if numEpisodeEdges>n
                        break;
                    end
                end
            end
        end
        %%
        function append_(obj, DataArray)
            % Append experience into ExperienceBuffer
            for ct=1:length(DataArray)
                % make sure observations and actions are stored as cell
                % arrays
                data = iCellifyExperience(DataArray{ct});
                % Saves data in memory
                obj.Memory{obj.NextIndex} = data;
                % Update isdone vector.
                obj.IsdoneVector(obj.NextIndex) = data{5};
                % Update isFull (useful later for drawing valid samples)
                if ~obj.IsFull && (obj.NextIndex == obj.Capacity)
                    obj.IsFull = true;
                end
                % Update NextIndex after data is stored
                if  obj.NextIndex == obj.Capacity
                    obj.NextIndex = 1;
                else
                    obj.NextIndex = obj.NextIndex + 1;
                end
            end
        end
        %%
        function [dataArray,idxes] = getSampledExperience_(obj,BatchSize,varargin)
            % Get sampled experiences for mini-batch to be created
            % dataArray = getSampledExperience(obj,BatchSize,gamma,n);
            if obj.Length < BatchSize
                dataArray = {};
                idxes = [];
            else
                idxes = randi(obj.Length, BatchSize, 1);
                dataArray = obj.Memory(idxes);
                % Using n-step look-ahead
                if nargin>2 && varargin{2}>1
                    dataArray = convertToNStepExperiences(obj,idxes,varargin{1},varargin{2});
                end
            end
        end
        %%
        function dataArray = getExperienceWithReturn_(obj,gamma,varargin)
            % Get experiences with return for mini-batch to be created
            % dataArray = getExperienceWithReturn(obj,gamma,n);
            
            [NewestDataIndex,OldestDataIndex] = findNewOldDataIdxes(obj);
            if obj.Memory{NewestDataIndex}{5}==0
                dataArray = {};
            else
                % Update Indexes for episodes to be learned
                n = 1;
                if nargin>2 && varargin{1}>1
                    n = varargin{1};
                end
                % REVISIT: use findEpisodeIsdoneIdx for better performance
                maxIdx = NewestDataIndex + (NewestDataIndex < OldestDataIndex) * obj.Capacity;
                [idxes,numEpisodeEdges] = findEpisodeIdxes(obj,maxIdx,OldestDataIndex,n);
                if numEpisodeEdges <= n
                    dataArray = {};
                    return;
                end
                dataArray = convertToMonteCarloExperiences(obj,idxes,gamma,n);
            end
        end
        %%
        function validateExperience(obj,NewExperience)
            % REVISIT: support batch experience
            % REVISIT: throw different error at each validation step
            
            % experience is a 1*5 cell array.
            [~,num] = size(NewExperience);
            if num~=5
                error((message('rl:general:errExperienceNumElements')));
            end
            Msg = message('rl:general:errExperienceSizeOrType');
            
            % validate dimension of NewExperience with data specs
            % validate 1st and 4th cell with ObsInfo
            NewObs = NewExperience{1};
            NewNextObs = NewExperience{4};
            
            assert(iIsDataTypeValid(NewObs), Msg);
            assert(iIsDataTypeValid(NewNextObs), Msg)
            for ct = 1:numel(obj.ObservationDimension)
                assert(iIsDataSizeValid(obj.ObservationDimension(ct), NewObs(ct)), Msg);
                assert(iIsDataSizeValid(obj.ObservationDimension(ct), NewNextObs(ct)), Msg);
            end
            % validate 2nd cell with ActInfo
            NewAction = NewExperience{2};
            assert(iIsDataTypeValid(NewAction), Msg);
            for ct = 1:numel(obj.ActionDimension)
                assert(iIsDataSizeValid(obj.ActionDimension(ct), NewAction(ct)), Msg);
            end
            % validate reward is scalar
            assert(isscalar(NewExperience{3}), Msg);
            % validate isdone is scalar, member of 0,1,2
            NewIsDone = NewExperience{5};
            assert(isscalar(NewIsDone) && ismember(NewIsDone, [0 1 2]), Msg);
        end
        %%
        function BatchExperience = batchExperience(obj, DataArray)
            % Concatenate 1xBatchSize of raw experience cell array into single
            % cell of Observation, Action, NextObservation, Reward, IsDone
            % Use observation and action specs dimension info for batching.
            
            % e.g: Each observation has dimension d x d x ... x d
            % {d x d x ... x d} x BatchSize => {d x d x ... x d x BatchSize}
            % Reward and IsDone: 1 x BatchSize
            
            if isempty(DataArray)
                BatchExperience = {};
            else
                BatchSize = numel(DataArray);
                
                % for-loop to address multiple input channels
                % e.g. image observation and vector observation
                for i = numel(DataArray{1}{1}):-1:1
                    % extract data from 1xN experience cell array
                    MiniBatchObservation{i}     = cellfun(@(x) x{1}{i},DataArray,'UniformOutput',false);
                    MiniBatchNextObservation{i} = cellfun(@(x) x{4}{i},DataArray,'UniformOutput',false);
                    sz = [obj.ObservationDimension{i} BatchSize];
                    
                    MiniBatchObservation{i}     = single(reshape([MiniBatchObservation{i}{:}],sz));
                    MiniBatchNextObservation{i} = single(reshape([MiniBatchNextObservation{i}{:}],sz));
                end
                for i = numel(DataArray{1}{2}):-1:1
                    % extract data from 1xN experience cell array
                    MiniBatchAction{i} = cellfun(@(x) x{2}{i},DataArray,'UniformOutput',false);
                    sz = [obj.ActionDimension{i} BatchSize];
                    
                    MiniBatchAction{i} = single(reshape([MiniBatchAction{i}{:}],sz));
                end
                MiniBatchReward = cellfun(@(x) single(x{3}),DataArray);
                MiniBatchIsDone = cellfun(@(x) single(x{5}),DataArray);
                
                BatchExperience = {MiniBatchObservation,MiniBatchAction,MiniBatchReward,MiniBatchNextObservation,MiniBatchIsDone};
            end
        end
        %%
        function [BatchExperience, Masks] = batchSequenceExperience(obj, DataArray, BatchSize, SequenceLength, Masks)
            % Concatenate 1x(SequenceLength*BatchSize) of raw experience 
            % cell array into single cell of Observation, Action, 
            % NextObservation, Reward, IsDone
            % Use observation and action specs dimension info for batching.
            
            % e.g: Each observation has dimension d x d x ... x d
            % {d x d x ... x d} x (SequenceLength*BatchSize) => {d x d x ... x d x BatchSize x SequenceLength}
            % Reward and IsDone: 1 x BatchSize x SequenceLength
            
            % If Masks is empty, it returns emtpy Masks.
            
            % Not support SequenceLenght = 1, use batchExperience instead
            
            if isempty(DataArray)
                BatchExperience = {};
                Masks = {};                                
            else                
                % for-loop to address multiple input channels
                % e.g. image observation and vector observation
                for i = numel(obj.ObservationDimension):-1:1
                    % extract data from 1x(SequenceLength*BatchSize) experience cell array
                    MiniBatchObservation{i}     = cellfun(@(x) x{1}{i},DataArray,'UniformOutput',false);
                    MiniBatchNextObservation{i} = cellfun(@(x) x{4}{i},DataArray,'UniformOutput',false);
                    DataSize = obj.ObservationDimension{i};
                    MiniBatchObservation{i} = iReshapeBatchSequence(...
                        MiniBatchObservation{i}, DataSize, BatchSize, SequenceLength);
                    MiniBatchNextObservation{i} = iReshapeBatchSequence(...
                        MiniBatchNextObservation{i}, DataSize, BatchSize, SequenceLength);
                end
                for i = numel(obj.ActionDimension):-1:1
                    % extract data from 1x(SequenceLength*BatchSize) experience cell array
                    MiniBatchAction{i}          = cellfun(@(x) x{2}{i},DataArray,'UniformOutput',false);
                    DataSize = obj.ActionDimension{i};
                    MiniBatchAction{i} = iReshapeBatchSequence(...
                        MiniBatchAction{i}, DataSize, BatchSize, SequenceLength);
                end
                
                % Reward + IsDone: 1xBxT
                DimsVec = [1 3 2];

                MiniBatchReward = cellfun(@(x) single(x{3}),DataArray);
                MiniBatchReward = reshape(MiniBatchReward,[1 SequenceLength BatchSize]);
                MiniBatchReward = permute(MiniBatchReward,DimsVec);

                MiniBatchIsDone = cellfun(@(x) single(x{5}),DataArray);
                MiniBatchIsDone = reshape(MiniBatchIsDone,[1 SequenceLength BatchSize]);
                MiniBatchIsDone = permute(MiniBatchIsDone,DimsVec);
                
                if ~isempty(Masks)
                    % Off-policy RL with a recurrent neural network uses
                    % masking to padd the inputs. 
                    Masks = reshape(Masks,[1 SequenceLength BatchSize]);
                    Masks = permute(Masks,DimsVec);
                end
                    
                BatchExperience = {MiniBatchObservation,MiniBatchAction,MiniBatchReward,MiniBatchNextObservation,MiniBatchIsDone};
            end
        end
    end
    
    methods (Static)
        function obj = loadobj(s)
            if isstruct(s)
                if ~isfield(s,'Version')
                    % version 1 does not have Version field
                    % In version 2,
                    %   - ValidateInputArguments name changes to DoValidate
                    %   - ObservationDimension and ActionDimension: cell 
                    %   arrays of observation and action dimensions
                    s.DoValidate = s.ValidateInputArguments;
                    s = rmfield(s,'ValidateInputArguments');
                    if isempty(s.Memory)
                        s.IsdoneVector = sparse(1, s.Capacity);
                        s.ObservationDimension = [];
                        s.ActionDimension      = [];
                    else
                        % collect ObservationDimension and ActionDimension
                        % from an experience in the memory
                        dummyExp = s.Memory{1};
                        [s.ObservationDimension, s.ActionDimension] = iGetDimInfoFromExp(dummyExp);
                        % reconstruct IsdoneVector
                        s.IsdoneVector = sparse(1, s.Capacity);
                        if s.IsFull
                            IsDoneVector = cellfun(@(x) single(x{5}), s.Memory);
                        else
                            IsDoneVector = cellfun(@(x) single(x{5}), s.Memory(1:s.NextIndex-1));
                        end
                        s.IsdoneVector(logical(IsDoneVector)) = IsDoneVector(logical(IsDoneVector));
                    end
                end
                
                obj = rl.util.ExperienceBuffer(s.Capacity,...
                    s.ObservationDimension,s.ActionDimension);
                obj.DoValidate = s.DoValidate;
                obj.IsFull = s.IsFull;
                obj.IsdoneVector = s.IsdoneVector;
                obj.Memory = s.Memory;
                obj.NextIndex = s.NextIndex;
                obj.SaveMemoryWithBuffer = s.SaveMemoryWithBuffer;
            else
                obj = s;
            end
            if isempty(obj.Memory)
                reset(obj);
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Local Functions%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function exp = iCellifyExperience(exp)
% observations and actions will be stored as a cell array
exp{1} = iCellify(exp{1});
exp{2} = iCellify(exp{2});
exp{4} = iCellify(exp{4});
end
function val = iCellify(val)
if ~iscell(val)
    val = {val};
end
end

function Data = iReshapeBatchSequence(Data, DataSize, BatchSize, SequenceLength)
% Not support SequenceLenght = 1

if iscell(Data)
    Data = single([Data{:}]);
end
% Data has size Cx(T*B)
if BatchSize == 1
    % reshape to Cx1xT
    Data = reshape(Data,[DataSize BatchSize SequenceLength]);
else % both BatchSize and SequenceLength > 1
    % reshape to CxTxB
    Data = reshape(Data,[DataSize SequenceLength BatchSize]);
    LastDim = ndims(Data);
    % swap last 2 dimension to get CxBxT
    DimsVec = [1:numel(DataSize) LastDim LastDim-1];
    % reshape to BxT
    Data = permute(Data,DimsVec);
end

end

function isValid = iIsDataSizeValid(DataDim,Data)
% checks if the Data has consistent size with DataDim specifications

isValid = false;
if ~iscell(Data)
    Data = {Data};
end
for ct = 1:numel(Data)
    isValid(ct) = isequal(DataDim{ct}, size(Data{ct}));
end
isValid = all(isValid);
end

function isValid = iIsDataTypeValid(Data)
% checks if the Data (Obs,NextObs,Act) has numeric data
% Data must be a cell array of numeric

isValid = true;
isValid = isValid && iscell(Data);
for ct = 1:numel(Data)
    isValid = isValid && isnumeric(Data{ct});
end
end

function [ObservationDimension, ActionDimension] = iGetDimInfoFromExp(SingleExp)
    % Extract observation and action dimension from a single experience
    % cell
    
    ObservationDimension = cellfun(@size,SingleExp{1},'UniformOutput',false);
    ActionDimension = cellfun(@size,SingleExp{2},'UniformOutput',false);
end