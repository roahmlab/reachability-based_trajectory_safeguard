function blkStruct = slblocks
% SLBLOCKS  Defines the Simulink library block representation
%           for Simulink Control Design

% Copyright 2017-2018 The MathWorks Inc.

blkStruct.Name    = sprintf('Reinforcement\nLearning');
blkStruct.OpenFcn = 'rllib';
blkStruct.MaskInitialization = '';
blkStruct.MaskDisplay = '';

% Define the library list for the Simulink Library browser.
% Return the name of the library model and the name for it.
Browser(1).Library = 'rllib';
Browser(1).Name    = 'Reinforcement Learning';
Browser(1).IsFlat  = 1; % Is this library "flat" (i.e. no subsystems)?

blkStruct.Browser = Browser;

% End of slblocks.m
