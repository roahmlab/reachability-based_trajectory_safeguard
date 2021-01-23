function rlCreateEnvTemplate(varargin)
% RLCREATEENVTEMPLATE: Create a template RL environment
%
%   RLCREATEENVTEMPLATE(CUSTOMCLASSNAME) creates a template class for RL 
%   environment and opens the associated CUSTOMCLASSNAME.m file. The 
%   template class contains a minimum implementation of a cart-pole 
%   environment. 
%   To define a custom environment, modify the template class.
%
%   Example
%       rlCreateEnvTemplate('myClass');
%
%   See also: rlFunctionEnv

%   Copyright 2017-2018 The MathWorks, Inc.

rl.util.rlCreateEnvTemplate(varargin{:});

end