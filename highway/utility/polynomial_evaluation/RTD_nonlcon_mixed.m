function [n, neq,gn,gneq] = RTD_nonlcon_mixed(x,wkcoef,wkpows,Jcoef,Jpows,N,timer,timeout)

    % get k without s
    keval = x;

    %ensure w(x,k)<1 for all points
    n=eval_w_4modes(keval,wkcoef,wkpows,N)-1;
    gn=eval_J_2Dk_4modes(keval,Jcoef,Jpows,N)';

    neq = [] ;

    gneq = [];
      if nargin>=7
       
        if toc(timer)>timeout
            
            error(['timeout',num2str(toc(timer))])
        end
    end
end


