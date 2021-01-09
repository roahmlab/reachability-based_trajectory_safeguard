function [n, neq,gn,gneq] = RTD_nonlcon_dynamic(x,wkcoef,wkpows,Jcoef,Jpows,N,timer,timeout)
    % get k without s
    keval = x;

    %ensure w(x,k)<1 for all points
    n=eval_w_Nmodes(keval,wkcoef,wkpows,N)-1;
    gn=eval_J_2Dk_Nmodes(keval,Jcoef,Jpows,N)';

    neq = [] ;
    gneq = [];
    
    if nargin>=7 
        if toc(timer)>timeout   
            error(['timeout',num2str(toc(timer))])
        end
    end
end


