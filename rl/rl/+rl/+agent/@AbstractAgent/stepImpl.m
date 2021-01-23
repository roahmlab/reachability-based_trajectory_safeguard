function action = stepImpl(this,exp)
% in training   mode: learn & return action with exploration
% in simulation mode: return action without exploration in simulation mode
% exp = {state,action,reward,nextstate,isdone}

% step agent, computing the appropriate action if training or not

% Copyright 2018 The MathWorks, Inc.

if ~isempty(this.CustomRewardFcn)
    exp{3} = feval(this.CustomRewardFcn,exp);
end

% evaluate q0
if this.EpisodeInfo.StepsTaken < 2
    this.EpisodeInfo.Q0 = evaluateQ0(this,exp);
end

% execute the function handle
action = this.StepFcn(this,exp);

end