function Learnables = scaleLearnables(Learnables, Scale)
% SCALELEARNABLES   Scale learnables. Learnables can be a cell array,
% numeric, or a table (Learnables property of dlnetwork).
%
%   Learnables = rl.internal.dataTransformation.scaleLearnables(Learnables, Scale) 
%   scales each learnables parameters in Learnables by Scale.

% Copyright 2019 The MathWorks Inc.

Learnables = dlupdate(@(x) x.*Scale, Learnables);
end