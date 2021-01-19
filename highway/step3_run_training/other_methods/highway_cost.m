function [cost, dcost]=highway_cost( K,kv,kh,x_des,start_tic,timeout)
% find derivative with respect to variable, here it is k peak and ky
% since ref position as a function of K has a closed form solution, the derivative also, given a desired x and y, this is for fmincon Naive RTD

%get final pos, subtrack initial pose, get second norm.
[x,y,~] = predict_highway_position( K,kv,kh );
% dy = alpha/120*5*tpk.^4+beta/24*4*tpk.^3+gam/6*3*tpk.^2 +a0/2*2*tpk+v0;
%turns out 
cost = sum((x_des -[x;y]).^2);
dcost = zeros(2,1);
%dc/dx* dx/dpk
dcost(1) = -2*(x_des(1)-x)*1;
dcost(2) = -2*(x_des(2)-y)*19/48;
if toc(start_tic) > timeout
        error('Timed out while evaluating cost function!')
end
%Also need cost gradient with respect to K, not x and y
