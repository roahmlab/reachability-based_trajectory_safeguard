function argStruct = generateDiscreteCriticMultiOutQFcn(argStruct,finiteSetSpec)
% GENERATEDISCRETECRITICMULTIOUTQFCN
% modifies argStruct to generate the body for critic based policies with
% multiple output Q critic

% Copyright 2019 The MathWorks, Inc.

% REVISIT: should this method live in Q representation instead of DQN agent

outputstr    = argStruct.OutputString;
inputstr     = argStruct.InputString ;
outputstr_y  = rl.codegen.handleMultipleOutputStrings(outputstr);

[elementstr,IL,IR] = rl.codegen.generateFiniteElementStrings(finiteSetSpec);
argStruct.EvaluateFunctionName = 'localEvaluate';

bodystr = sprintf([...
    'actionSet = %s;\n',...
    'q = localEvaluate(%s);\n',...
    '[~,actionIndex] = max(q);\n',...
    '%s = actionSet%sactionIndex%s;'],...
    elementstr,inputstr,outputstr_y,IL,IR...
    );
argStruct.EvaluateOutputString = 'q';
argStruct.EvaluateInputString = inputstr;
argStruct.BodyString = bodystr;

