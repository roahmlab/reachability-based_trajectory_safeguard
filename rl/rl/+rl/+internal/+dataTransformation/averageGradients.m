function gavg = averageGradients(gs)
% AVERAGEGRADIENTS Averages a set of gradients in structure format (e.g.
% s.Critic, s.Actor).

% Copyright 2019 The MathWorks Inc

% assumes the grads are in structure format (e.g. s.Critic, s.Actor)

% get the number of gradients to average
N = numel(gs);
fs = fields(gs);
gavg = gs(1);
for i = 1:numel(fs)
    % loop through each field
    f = fs{i};
    % elementwise sum all gradients
    for j = 2:N
        gavg.(f) = gadd(gavg.(f),gs(j).(f));
    end
    % elementwise divide all gradient by N
    gavg.(f) = gdivide(gavg.(f),N);
end