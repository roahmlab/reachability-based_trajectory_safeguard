function Options = rlDQNAgentOptions(varargin)
% rlDQNAgentOptions: Creates options for DQN Agent.
%
%   OPT = rlDQNAgentOptions returns the default options for rlDQNAgent. 
%
%   OPT = rlDQNAgentOptions('Option1',Value1,'Option2',Value2,...) uses name/value
%   pairs to override the default values for 'Option1','Option2',...
%
%   Supported options are:
%
%   UseDoubleDQN                        Use Double DQN for learning
%   EpsilonGreedyExploration            Parameters for Epsilon Greedy exploration
%       Epsilon                         Probability threshold for agent to either randomly
%                                       select a valid action or select the action that 
%                                       maximizes the state-action value function
%       EpsilonMin                      Minimum value of Epsilon
%       EpsilonDecay                    Decay rate of Epsilon when Epsilon is updated
%   MiniBatchSize                       Size of the mini-batch of experiences sampled from the experience 
%                                       buffer for each learning iteration.
%                                       If the critic network is a
%                                       recurrent neural network (LSTM for
%                                       example), ''MiniBatchSize'' is the
%                                       number of the trajectories in a batch.
%   TargetSmoothFactor                  Smoothing factor that determines how a target model used for
%                                       prediction is updated with weights from the prediction model used
%                                       for training
%   TargetUpdateFrequency               Number of steps after which the target model is updated
%   ResetExperienceBufferBeforeTraining Empty the experience buffer when training begins
%   NumStepsToLookAhead                 Number of steps to look-ahead when computing return
%                                       N-step Q learning is not supported for a recurrent neural network.
%                                       If a reccurent neural network is used, NumStepsToLookAhead must be 1.
%   ExperienceBufferLength              The length of the buffer to store experiences used to
%                                       sample mini-batch of experiences in training
%   SampleTime                          Sample time of the agent
%   DiscountFactor                      Discount factor to apply to future rewards during training
%   SaveExperienceBufferWithAgent       Save the experience buffer when saving the agent
%   SequenceLength                      Maximum trajectory length in batch training when a recurrent neural
%                                       network is used for the critic network
%
%   See also: rlDQNAgent, rlDDPGAgentOptions, rlPGAgentOptions, rlACAgentOptions

% Copyright 2017-2018 The MathWorks Inc.

Options = rl.option.rlDQNAgentOptions(varargin{:});

end