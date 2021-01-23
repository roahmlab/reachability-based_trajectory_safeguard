function trainer = createTrainerFactory(env,agent,trainingOptions)
% CREATETRAINERFACTORY
%
% Create a trainer object given env, agent, and training options

% Revised: 8-22-2019
% Copyright 2019 The MathWorks, Inc.

if trainingOptions.UseParallel
    
    popts           = trainingOptions.ParallelizationOptions;
    parallelization = popts.Mode;
    data2send       = popts.DataToSendFromWorkers;
    sendexp         = strcmpi(data2send,"experiences");
    
    switch parallelization
        case 'async'
            if sendexp
                trainer = rl.train.dq.ExpParCommTrainer(env,agent,trainingOptions);
            else
                trainer = rl.train.dq.GradParCommTrainer(env,agent,trainingOptions);
            end
        case 'sync'
            % use parfor if we have to wait till the end of each episode.
            % Otherwise use the comm interface.
            stepsUntilDataIsSent = popts.StepsUntilDataIsSent;
            if stepsUntilDataIsSent == -1
                trainer = rl.train.parfor.ParforTrainer(env,agent,trainingOptions);
            else
                if sendexp
                    trainer = rl.train.dq.ExpParCommTrainer(env,agent,trainingOptions);
                else
                    trainer = rl.train.dq.GradParCommTrainer(env,agent,trainingOptions);
                end
            end
        otherwise
            % this should never be hit
            error(message('rl:general:TrainingManagerInvalidTrainingScheme',parallelization));
    end
else
    trainer = rl.train.SeriesTrainer(env,agent,trainingOptions);
end

