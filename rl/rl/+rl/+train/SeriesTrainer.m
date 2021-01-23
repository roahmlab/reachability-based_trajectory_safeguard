classdef SeriesTrainer < rl.train.Trainer
% SERIESTRAINER
%
% Train an agent against an environment on a single thread

% Revised: 7-2-2019
% Copyright 2019 The MathWorks, Inc.

    properties (Access = private)
        Listeners = event.listener.empty
    end
    methods
        function this = SeriesTrainer(env,agent,trainOpts)
            this = this@rl.train.Trainer(env,agent,trainOpts);
            
            this.Listeners(1) = addlistener(env,'EpisodeFinished',...
                    @(src,ed) notifyEpisodeFinishedAndCheckStopTrain(this,ed.Data));
        end
        function delete(this)
            delete(this.Listeners)
        end
        function run(this)
            seriestaskspec = rl.task.SeriesTrainTaskSpec(this.Env,this.Agent,this.TrainOpts);
            run(seriestaskspec);
        end
    end
end