%% description
% This script creates a quadrotor_agentHelper class instance and uses its
% adjust method to make sure a randomly-generated trajectory is safe. It
% then plots the agentHelper.
%
% Author: Shreyas Kousik
% Created: 2 Aug 2020
% Updated: 4 Aug 2020

%% create agentHelper
% get agent, FRS, and tracking error table
A = quadrotor_agent('camera_view_style','behind') ;
FRS = load('quadrotor_FRS_v7_a10_dt0.02.mat') ;
if ~exist('tbl','var')
    tbl = load('quadrotor_tracking_error_table_dt0.02_vmax_5.25_zonotope.mat') ;
else
    disp('table already loaded') ;
end

% create agentHelper
AH = quadrotor_agentHelper(A,FRS,tbl,'v_max',5,'a_max',3,'verbose',10) ;

% create world
W = zonotope_box_world('start',[2;0;2]) ;

% reset agent
AH.reset([],W.start) ;

%% test if a random k_user is unsafe
% create random k_user
k_user = 10*rand(3,1) - 5 ;

% test if it is unafe
AI = AH.A.get_agent_info() ;
WI = W.get_world_info(AI) ;
k = AH.adjust(k_user,WI) ;

%% plotting
figure(1) ; clf ; axis equal ; hold on ;

plot(W) ;
plot(AH) ;