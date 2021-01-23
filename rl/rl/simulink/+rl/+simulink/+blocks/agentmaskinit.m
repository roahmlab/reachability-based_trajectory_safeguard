function agentmaskinit(blk)
% AGENTMASKINIT: Initialize RL agent mask.

% Copyright 2018 The MathWorks, Inc.

% get the root level model
rootmdl = bdroot(blk);

% set the wrapper block properties
wrapperblk = [blk,'/AgentWrapper'];
set_param(wrapperblk,'ModelName',rootmdl   );
set_param(wrapperblk,'BlockPath',wrapperblk);

set_param(wrapperblk,'IsRewardProvided',localOnOff2TrueFalse(get_param(blk,'ProvideReward')));
set_param(wrapperblk,'IsIsDoneProvided',localOnOff2TrueFalse(get_param(blk,'ProvideIsDone')));

set_param(wrapperblk,'CollectLoggedSignals',localOnOff2TrueFalse(get_param(blk,'CollectLoggedSignals')));
set_param(wrapperblk,'CollectModelStates'  ,localOnOff2TrueFalse(get_param(blk,'CollectModelStates'  )));

set_param(wrapperblk,'TreatAsDirectFeedthrough',localOnOff2TrueFalse(get_param(blk,'TreatAsDirectFeedthrough')));

% update the ports based on selections
rl.simulink.blocks.agentmaskportupdate(blk);

%% Local functions
function val = localOnOff2TrueFalse(val)
if strcmp(val,'on')
    val = 'true';
else
    val = 'false';
end

