function opt = rlSACAgentOptions(varargin)
%rlSACAgentOptions: Creates options for SAC Agent
%
%   OPT = rlSACAgentoptions returns the default options for rlSACAgent. 
%
%   OPT = rlSACAgentOptions('Option1',Value1,'Option2',Value2,...) uses name/value
%   pairs to override the default values for 'Option1','Option2',...
%
%   Supported options are:
%
%   EntropyWeightOptions                Options for tuning Entropy Weight
%   PolicyUpdateFrequency               Number of steps after which the policy is updated
%   TargetUpdateFrequency               Number of steps after which the target model is updated
%   TargetSmoothFactor                  Smoothing factor that determines how a target model used for
%                                       prediction is updated with weights from the prediction model used
%                                       for training
%   MiniBatchSize                       Size of the mini-batch of experiences sampled from the experience 
%                                       buffer for each learning iteration
%   ResetExperienceBufferBeforeTraining Empty the experience buffer when training begins
%   NumStepsToLookAhead                 Number of steps to look-ahead when computing return
%   ExperienceBufferLength              The length of the buffer to store experiences used to
%                                       sample mini-batch of experiences in training
%   SampleTime                          Sample time of the agent
%   DiscountFactor                      Discount factor to apply to future rewards during training
%   SaveExperienceBufferWithAgent       Save the experience buffer when saving the agent
%
%   See also: rlSACAgent, EntropyWeightOptions, rlTD3AgentOptions, rlPPOAgentOptions, rlDDPGAgentOptions

% Copyright 2019 The MathWorks, Inc.

opt = rl.option.rlSACAgentOptions(varargin{:});
end