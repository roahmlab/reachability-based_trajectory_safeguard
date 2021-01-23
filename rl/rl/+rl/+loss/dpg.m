function Loss = dpg(OutputAction,LossVariable)
% Deterministic policy loss.
%   OutputAction: action output from deterministic actor representation
%   LossVariable: struct contains observations and the agent's Critic

% Copyright 2019 The MathWorks Inc.

ObsDim = finddim(OutputAction,'B');
NumObs = size(OutputAction, ObsDim);
for ct = 1:numel(LossVariable.Critic)
    if ct < 2
        CriticOutput = getValue(LossVariable.Critic(ct), LossVariable.Observation, {OutputAction});
    else
        CriticOutput = min(CriticOutput,getValue(LossVariable.Critic(ct), LossVariable.Observation, {OutputAction}));
    end
end
Loss = -sum(CriticOutput)/NumObs;
end