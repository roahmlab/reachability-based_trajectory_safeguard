function guihelp(topickey,varargin)
% GUIHELP

% Revised: 10-30-2018
% Copyright 2018 The MathWorks, Inc.

helpdir = docroot;
if ~isempty(helpdir)
    mapfile = fullfile(helpdir,'reinforcement-learning','helptargets.map');
else
    mapfile = '';
end
helpview(mapfile,topickey);

