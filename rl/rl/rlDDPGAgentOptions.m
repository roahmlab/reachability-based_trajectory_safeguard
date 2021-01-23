function opt = rlDDPGAgentOptions(varargin)
%rlDDPGAgentOptions: Creates options for DDPG Agent
%
%   OPT = rlDDPGAgentOptions returns the default options for rlDDPGAgent. 
%
%   OPT = rlDDPGAgentOptions('Option1',Value1,'Option2',Value2,...) uses name/value
%   pairs to override the default values for 'Option1','Option2',...
%
%   Supported options are:
%
%   NoiseOptions                        Parameters for Ornstein Uhlenbeck noise
%       InitialAction                       Initial state of the noise model
%       Mean                                Mean of the noise model
%       MeanAttractionConstant              Constant used to attract the process toward the mean
%       Variance                            Variance of the random process
%       VarianceDecayRate                   Rate of noise variance decay with each step of the noise model 
%   MiniBatchSize                       Size of the mini-batch of experiences sampled from the experience 
%                                       buffer for each learning iteration
%   TargetSmoothFactor                  Smoothing factor that determines how a target model used for
%                                       prediction is updated with weights from the prediction model used
%                                       for training
%   TargetUpdateFrequency               Number of steps after which the target model is updated
%   ResetExperienceBufferBeforeTraining Empty the experience buffer when training begins
%   NumStepsToLookAhead                 Number of steps to look-ahead when computing return
%   ExperienceBufferLength              The length of the buffer to store experiences used to
%                                       sample mini-batch of experiences in training
%   SampleTime                          Sample time of the agent
%   DiscountFactor                      Discount factor to apply to future rewards during training
%   SaveExperienceBufferWithAgent       Save the experience buffer when saving the agent
%
%   See also: rlDDPGAgent, rlDQNAgentOptions, rlPGAgentOptions, rlACAgentOptions

% Copyright 2018 The MathWorks, Inc.

opt = rl.option.rlDDPGAgentOptions(varargin{:});
end
