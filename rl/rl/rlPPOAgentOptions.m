function Options = rlPPOAgentOptions(varargin)
% rlPPOAgentOptions: Creates options for PPO Agent.
%
%   OPT = rlPPOAgentOptions returns the default options for rlPPOAgent. 
%
%   OPT = rlPPOAgentOptions('Option1',Value1,'Option2',Value2,...) uses name/value
%   pairs to override the default values for 'Option1','Option2',...
%
%   Supported options are:
%
%   ExperienceHorizon               Number of steps an agent interacts with the environment before it learns from its experience.
%   ClipFactor                      Clip factor to limit the change in each policy update step.
%   EntropyLossWeight               Weight for entropy loss to promote policy exploration.
%   MiniBatchSize                   The size of the mini-batch used for each training iteration.
%   NumEpoch                        Number of epochs the actor and critic learn from the current experience set.
%   AdvantageEstimateMethod         Method used for estimating advantage values ('finite-horizon' or 'gae').
%   GAEFactor                       Smoothing factor of the generalized advantage estimator.
%   SampleTime                      Sample time of the agent.
%   DiscountFactor                  Discount factor to apply to future rewards during training.
%
%   See also: rlPPOAgent, rlDDPGAgentOptions, rlPPOAgentOptions, rlDQNAgentOptions

% Copyright 2019 The MathWorks, Inc.

Options = rl.option.rlPPOAgentOptions(varargin{:});

end