classdef rlTrainingOptions < rl.option.BaseSimOptions
    % RLTRAININGOPTIONS: Creates training options for RL.
    %
    %   OPT = RLTRAININGOPTIONS returns the default options for RL agents.
    %
    %   OPT = RLTRAININGOPTIONS('Option1',Value1,'Option2',Value2,...) uses name/value
    %   pairs to override the default values for 'Option1','Option2',...
    %
    %   Supported options are:
    %
    %   MaxEpisodes                 Maximum number of episodes to train the agent
    %   MaxStepsPerEpisode          Maximum number of steps to run per episode
    %   ScoreAveragingWindowLength  Window length for averaging scores and rewards
    %   StopTrainingCriteria        Stop criteria for agent training
    %                               - AverageSteps   : Running average of number of steps per episode
    %                               - AverageReward  : Running average of reward per episode
    %                               - EpisodeReward  : Reward for current episode
    %                               - GlobalStepCount: Total times the agent was invoked
    %                               - EpisodeCount   : Total number of episodes the agent has trained for
    %   StopTrainingValue           Value of training criteria to stop
    %   SaveAgentCriteria           Criteria to save agent while training
    %                               - none           : Do not save any agent in training
    %                               - AverageSteps   : Running average of number of steps per episode
    %                               - AverageReward  : Running average of reward per episode
    %                               - EpisodeReward  : Reward for current episode
    %                               - GlobalStepCount: Total times the agent was invoked
    %                               - EpisodeCount   : Total number of episodes the agent has trained for
    %   SaveAgentValue              Value of criteria to save
    %   SaveAgentDirectory          Directory name for saving agents
    %   UseParallel                 Use parallel training schemes defined by ParallelizationOptions (true, false (default)
    %   ParallelizationOptions      Parallelization options to control parallel training schemes
    %                               - Mode                           : Parallel training scheme
    %                                   - sync (default)             : Use the parpool to run synchronous training on the available workers
    %                                   - async                      : Use the parpool to run asynchronous training on the available workers
    %                               - DataToSendFromWorkers          : Data that will be sent from the workers to the client ("experiences" (default), "gradients").
    %                               - StepsUntilDataIsSent           : The number of steps the worker will take before sending data back to the client. If -1 (default), the worker will send data at the end of an episode.
    %                               - WorkerRandomSeeds              : Initialize the random seed on the worker using <a href="matlab: help rng">rng</a>. If -1 (default) the random seed on each worker will be initialized with the worker ID. Use a vector to initialize a seed for each worker.
    %                               - TransferBaseWorkspaceVariables : Transfer base workspace variables to each worker before training starts ("on" (default), "off").
    %                               - AttachedFiles                  : Files to attach to the parpool
    %                               - SetupFcn                       : Setup function handle to run on each worker before training starts.
    %                               - CleanupFcn                     : Cleanup function handle to run on each worker after training ends.
    %   StopOnError                 Stop remaining simulations if any simulation errors occurs. ("on" (default), "off")
    %   Verbose                     Display training progress on the command line (false (default), true)
    %   Plots                       Display training progress with the Episode Manager ("training-progress" (default), "none")
    %
    % See also: <a href="matlab: help rl.agent.AbstractAgent\train">train</a>, parpool

    % Copyright 2017-2018 The MathWorks Inc.    
    
    properties (Access = protected)
        StopTrainingCriteria_
        StopTrainingValue_
        
        SaveAgentCriteria_
        SaveAgentValue_
    end
    
    properties
        % Maximum number of episodes to train the agent
        MaxEpisodes
        
        % Maximum number of steps to run per episode
        MaxStepsPerEpisode
        
        % Window length for averaging scores and rewards
        ScoreAveragingWindowLength                
    end
    
    properties (Dependent)
        % Stop criteria for training and saving
        % Possible options
        % - AverageSteps   : Running average of number of steps per episode
        % - AverageReward  : Running average of reward per episode
        % - EpisodeReward  : Reward for current episode
        % - GlobalStepCount: Total times the agent was invoked
        % - EpisodeCount   : Total number of episodes the agent has trained for
        StopTrainingCriteria
        StopTrainingValue
        
        SaveAgentCriteria
        SaveAgentValue   
    end
    
    properties
        % directory to save agents
        SaveAgentDirectory
        
        % Print detailed information in command window
        Verbose  
        
        % Plots to display training progress
        Plots
    end  
    properties (Dependent,Hidden)
        % for existing compatibility
        Parallelization
    end
    
    methods
        function obj = rlTrainingOptions(varargin)
            %% Parsing name-value pairs
            
            obj = obj@rl.option.BaseSimOptions(varargin{:});
            
            parser = inputParser;
            obj.addParseParams(parser);                        
            parse(parser,varargin{:});             
            
            %% Property assignments                                
            obj.MaxEpisodes = parser.Results.MaxEpisodes;
            obj.MaxStepsPerEpisode = parser.Results.MaxStepsPerEpisode;
            obj.ScoreAveragingWindowLength = parser.Results.ScoreAveragingWindowLength;
            
            obj.Verbose = parser.Results.Verbose;                     
            obj.Plots = parser.Results.Plots;                     
                        
            if isempty(parser.Results.StopTrainingCriteria) && isa(parser.Results.StopTrainingCriteria,'double')
                % no need to validate if default value is set
                obj.StopTrainingCriteria_ = "AverageSteps";
            else
                % validate stop training criteria when user provides
                obj.StopTrainingCriteria_ = rl.option.rlTrainingOptions.validateStopTrainingCriteria(parser.Results.StopTrainingCriteria);
            end
            
            if isempty(parser.Results.StopTrainingValue)
                % no need to validate if default value is set
                obj.StopTrainingValue_ = rl.option.rlTrainingOptions.defaultCriteriaValue(obj.StopTrainingCriteria);
            else
                % validate stop training value when user provides
                obj.StopTrainingValue = parser.Results.StopTrainingValue;
            end
            
            rl.option.rlTrainingOptions.validateCriteriaValue(obj.StopTrainingCriteria,'StopTrainingCriteria',obj.StopTrainingValue)
            
            if isempty(parser.Results.SaveAgentCriteria) && isa(parser.Results.SaveAgentCriteria,'double')
                % no need to validate if default value is set
                obj.SaveAgentCriteria_ = "none";
            else
                % validate save agent criteria when user provides
                obj.SaveAgentCriteria_ = rl.option.rlTrainingOptions.validateSaveAgentCriteria(parser.Results.SaveAgentCriteria);
            end            
                        
            if isempty(parser.Results.SaveAgentValue)
                % no need to validate if default value is set
                obj.SaveAgentValue_ = rl.option.rlTrainingOptions.defaultCriteriaValue(obj.SaveAgentCriteria);
            else
                % validate save agent value when user provides
                obj.SaveAgentValue = parser.Results.SaveAgentValue;
            end
            
            rl.option.rlTrainingOptions.validateCriteriaValue(obj.SaveAgentCriteria,'SaveAgentCriteria',obj.SaveAgentValue)
            
            obj.SaveAgentDirectory = parser.Results.SaveAgentDirectory;
        end
        function obj = set.MaxEpisodes(obj,Value)
            validateattributes(Value,{'numeric'},{'scalar','real','integer','positive'},'','MaxEpisodes');
            obj.MaxEpisodes = Value;
        end
        function obj = set.MaxStepsPerEpisode(obj,Value)
            validateattributes(Value,{'numeric'},{'scalar','real','integer','positive'},'','MaxStepsPerEpisode');
            obj.MaxStepsPerEpisode = Value;
        end
        function obj = set.ScoreAveragingWindowLength(obj,Value)
            validateattributes(Value,{'numeric'},{'scalar','real','integer','positive'},'','ScoreAveragingWindowLength');
            obj.ScoreAveragingWindowLength = Value;
        end
        function obj = set.Verbose(obj,Value)
            validateattributes(Value,{'logical','numeric'},{'scalar','real'},'','Verbose');
            obj.Verbose = logical(Value);
        end        
        function obj = set.Plots(obj,Value)
            Criteria = rl.option.rlTrainingOptions.getPlots;
            validatestring(Value,Criteria,'','Plots');
            obj.Plots = string(Value);
        end
        function obj = set.StopTrainingCriteria(obj,CriteriaValue)                  
            obj.StopTrainingCriteria_ = rl.option.rlTrainingOptions.validateStopTrainingCriteria(CriteriaValue);
            try
                rl.option.rlTrainingOptions.validateCriteriaValue(obj.StopTrainingCriteria,'StopTrainingValue',obj.StopTrainingValue);
            catch
                obj.StopTrainingValue = rl.option.rlTrainingOptions.defaultCriteriaValue(obj.StopTrainingCriteria);
            end
        end
        function criteria = get.StopTrainingCriteria(obj)
            criteria = obj.StopTrainingCriteria_;            
        end
        function obj = set.StopTrainingValue(obj,Value)
            rl.option.rlTrainingOptions.validateCriteriaValue(obj.StopTrainingCriteria,'StopTrainingValue',Value);
            obj.StopTrainingValue_ = Value;
        end
        function val = get.StopTrainingValue(obj)
            val = obj.StopTrainingValue_;            
        end        
        function obj = set.SaveAgentCriteria(obj,CriteriaValue)                   
            obj.SaveAgentCriteria_ = rl.option.rlTrainingOptions.validateSaveAgentCriteria(CriteriaValue);
            try
                rl.option.rlTrainingOptions.validateCriteriaValue(obj.SaveAgentCriteria,'SaveAgentValue',obj.SaveAgentValue);
            catch
                obj.SaveAgentValue = rl.option.rlTrainingOptions.defaultCriteriaValue(obj.SaveAgentCriteria);
            end            
        end
        function criteria = get.SaveAgentCriteria(obj)
            criteria = obj.SaveAgentCriteria_;            
        end        
        function obj = set.SaveAgentValue(obj,Value)
            rl.option.rlTrainingOptions.validateCriteriaValue(obj.SaveAgentCriteria,'SaveAgentValue',Value);
            obj.SaveAgentValue_ = Value;
        end
        function val = get.SaveAgentValue(obj)
            val = obj.SaveAgentValue_;            
        end         
        function obj = set.SaveAgentDirectory(obj,Value)
            validateattributes(Value,{'char','string'},{},'','SaveAgentDirectory');
            obj.SaveAgentDirectory = string(Value);
        end   
        function obj = set.Parallelization(obj,val)
            % for compatibility
            strval = validatestring(val,[{'none'},rl.option.ParallelTraining.ValidParallelModes],'','Parallelization');
            if strcmp(strval,'none')
                obj.UseParallel_ = false;
            else
                obj.UseParallel_ = true;
                obj.ParallelizationOptions_.Mode = strval;
            end
        end
        function val = get.Parallelization(obj)
            % for compatibility
            if obj.UseParallel_
                val = obj.ParallelizationOptions_.Mode;
            else
                val = "none";
            end
        end
        function opts = getSimulationOptions(obj)
            % cast the paropts
            paropts = rl.option.ParallelSimulation(obj.ParallelizationOptions);
            opts = rlSimulationOptions('MaxSteps',obj.MaxStepsPerEpisode,...
                'NumSimulations',obj.MaxEpisodes,...
                'StopOnError',obj.StopOnError,...
                'UseParallel',obj.UseParallel,...
                'ParallelizationOptions',paropts);
        end
    end
    methods (Hidden)        
        function flag = StopTrainingFunction(obj,info)
            if any(strcmpi(obj.StopTrainingCriteria,rl.option.rlTrainingOptions.getFunctionHandleCriteria))
                flag = obj.StopTrainingValue(info);
            else
                flag = info.(obj.StopTrainingCriteria)>=obj.StopTrainingValue;
            end
        end
        function flag = SaveAgentFunction(obj,info)
            if strcmpi(obj.SaveAgentCriteria,"none")
                flag = false;
            else
                if any(strcmpi(obj.SaveAgentCriteria,rl.option.rlTrainingOptions.getFunctionHandleCriteria))
                    flag = obj.SaveAgentValue(info);
                else
                    flag = info.(obj.SaveAgentCriteria)>=obj.SaveAgentValue;
                end
            end
        end
        function validateStopTrainingFunction(obj)
            criteria = rl.option.rlTrainingOptions.getCustomCriteria;
            info = cell2struct(num2cell(ones(1,numel(criteria))),criteria,2);
            try
                flag = StopTrainingFunction(obj,info);
            catch
                error(message('rl:general:TrainingOptionsErrorStopTrainingCustomFunction'));                
            end
            if ~( (isnumeric(flag) || islogical(flag)) && isscalar(flag) )
                error(message('rl:general:TrainingOptionsErrorStopTrainingCustomFunction'));
            end            
        end
        function validateSaveAgentFunction(obj)
            criteria = rl.option.rlTrainingOptions.getCustomCriteria;
            info = cell2struct(num2cell(ones(1,numel(criteria))),criteria,2);
            try
                flag = SaveAgentFunction(obj,info);
            catch
                error(message('rl:general:TrainingOptionsErrorSaveAgentCustomFunction'));
            end
            if ~( (isnumeric(flag) || islogical(flag)) && isscalar(flag) )
                error(message('rl:general:TrainingOptionsErrorSaveAgentCustomFunction'));
            end            
        end        
    end
    methods (Static, Access = private)
        function validateCriteriaValue(Criteria,CriteriaValueName,Value)
            switch Criteria
                case rl.option.rlTrainingOptions.getPositiveIntegerCriteria
                    validateattributes(Value,{'numeric'},{'scalar','real','integer','positive'},'',CriteriaValueName);
                case rl.option.rlTrainingOptions.getRealNumberCriteria
                    validateattributes(Value,{'numeric'},{'scalar','real','nonnan'},'',CriteriaValueName);  
                case rl.option.rlTrainingOptions.getFunctionHandleCriteria
                    if isa(Value,'function_handle') && exist(func2str(Value),'file')~=2
                        warning(message('rl:general:TrainingOptionsWarningCustomCriteria'))
                    end                    
                    validateattributes(Value,{'function_handle'},{'scalar'},'',CriteriaValueName);
                case rl.option.rlTrainingOptions.getnoneCriteria
                    if ~strcmpi(Value,"none")
                        error(message('rl:general:TrainingOptionsErrorSaveAgentValue'));
                    end
            end
        end
        function CurrentCriteria = validateStopTrainingCriteria(CriteriaValue)
            Criteria = rl.option.rlTrainingOptions.getStopTrainingCriteria;
            CurrentCriteria = string(validatestring(CriteriaValue,Criteria,'','StopTrainingCriteria'));
        end
        function CurrentCriteria = validateSaveAgentCriteria(CriteriaValue)
            Criteria = rl.option.rlTrainingOptions.getSaveAgentCriteria;
            CurrentCriteria = string(validatestring(CriteriaValue,Criteria,'','SaveAgentCriteria'));
        end            
        function Criteria = getPositiveIntegerCriteria()
            Criteria = {'AverageSteps', ...
                        'GlobalStepCount', ...
                        'EpisodeCount', ...
                        'EpisodeSteps'};
        end
        function Criteria = getRealNumberCriteria()
            Criteria = {'AverageReward', ...
                        'EpisodeReward'};
        end
        function Criteria = getFunctionHandleCriteria()
            Criteria = {'Custom'};
        end
        function Criteria = getnoneCriteria()
            Criteria = {'none'};
        end        
        function Criteria = getStopTrainingCriteria()
            Criteria = [ rl.option.rlTrainingOptions.getPositiveIntegerCriteria, ...
                         rl.option.rlTrainingOptions.getRealNumberCriteria, ...
                         rl.option.rlTrainingOptions.getFunctionHandleCriteria ];
        end
        function Criteria = getSaveAgentCriteria()
            Criteria = [ rl.option.rlTrainingOptions.getPositiveIntegerCriteria, ...
                         rl.option.rlTrainingOptions.getRealNumberCriteria, ...
                         rl.option.rlTrainingOptions.getFunctionHandleCriteria, ...
                         rl.option.rlTrainingOptions.getnoneCriteria];
        end                
        function Criteria = getCustomCriteria()
            Criteria = [ rl.option.rlTrainingOptions.getPositiveIntegerCriteria, ...
                         rl.option.rlTrainingOptions.getRealNumberCriteria];
        end                
        function Value = defaultCriteriaValue(Criteria)
            switch Criteria
                case rl.option.rlTrainingOptions.getPositiveIntegerCriteria
                    Value = 500;
                case rl.option.rlTrainingOptions.getRealNumberCriteria
                    Value = 500;
                case rl.option.rlTrainingOptions.getFunctionHandleCriteria
                    Value = @valueFunction;
                case rl.option.rlTrainingOptions.getnoneCriteria
                    Value = "none";
            end          
        end 
        function Plots = getPlots()
            Plots = {'none', 'training-progress'};
        end
    end
    methods (Access = protected,Static)
        function val = validateParallelizationOptions(val)
            validateattributes(val,{'rl.option.ParallelTraining'},{'scalar'},'','ParallelizationOptions');
        end
        function opts = getDefaultParallelOptions()
            opts = rl.option.ParallelTraining();
        end
    end
    methods (Access = protected)
        function addParseParams(obj,parser) 
            addParseParams@rl.option.BaseSimOptions(obj,parser);
            % Parameters for training the model
            addParameter(parser,'MaxEpisodes',500);
            addParameter(parser,'MaxStepsPerEpisode',500);
            addParameter(parser,'ScoreAveragingWindowLength',5);
                        
            % Parameter for display
            addParameter(parser,'Verbose',false);            
            addParameter(parser,'Plots','training-progress');            
            
            % Parameters for stop training criteria and threshold
            addParameter(parser,'StopTrainingCriteria',[]);
            addParameter(parser,'StopTrainingValue',[]);
            
            % Parameters for save agent criteria, threshold, and directory
            addParameter(parser,'SaveAgentCriteria',[]);
            addParameter(parser,'SaveAgentValue',[]);
            addParameter(parser,'SaveAgentDirectory',"savedAgents");
        end
    end
end