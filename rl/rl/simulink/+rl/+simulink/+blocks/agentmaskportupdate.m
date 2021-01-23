function agentmaskportupdate(blk)
% AGENTMASKPORTUPDATE

% Copyright 2018 The MathWorks, Inc.

prwd = strcmp(get_param(blk,'ProvideReward'),'on');
pisd = strcmp(get_param(blk,'ProvideIsDone'),'on');
pcms = strcmp(get_param(blk,'ProvideCumRwd'),'on');

path = getfullname(blk);
rwdblk = [path,'/reward'];
isdblk = [path,'/isdone'];
cmsblk = [path,'/cumulative_reward'];

% replace input port with ground or visa versa
localReplaceWithGround(rwdblk,prwd);
localReplaceWithGround(isdblk,pisd);

% replace output port with terminator or visa versa
localReplaceWithTerminator(cmsblk,pcms);

% make sure the reward port is port 2
if strcmp(get_param(rwdblk,'BlockType'),'Inport')
    set_param(rwdblk,'Port','2');
end

function localReplaceWithGround(innerblkpath,wantsport)
% replace the port with ground or input port
if wantsport
    newblk = 'simulink/Sources/In1';
    replace = ~strcmp(get_param(innerblkpath,'BlockType'),'Inport');
else
    newblk = 'simulink/Sources/Ground';
    replace = ~strcmp(get_param(innerblkpath,'BlockType'),'Ground');
end
if replace
    slInternal('replace_block',innerblkpath,newblk);
end

function localReplaceWithTerminator(innerblkpath,wantsport)
% replace the port with ground or input port
if wantsport
    newblk = 'simulink/Sinks/Out1';
    replace = ~strcmp(get_param(innerblkpath,'BlockType'),'Outport');
else
    newblk = 'simulink/Sinks/Terminator';
    replace = ~strcmp(get_param(innerblkpath,'BlockType'),'Terminator');
end
if replace
    slInternal('replace_block',innerblkpath,newblk);
end




