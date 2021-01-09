function [n, neq,gn,gneq] = RTD_nonlcon_nDK_nModes(x,wkcoef,wkpows,Jcoef,Jpows,N,timer,timeout)
    % get k without s
    keval = x;
    
    neq = [] ;
    gneq = [];
    
    if ~isempty(N)
        
        %ensure w(x,k)<1 for all points
        n=eval_w_Nmodes(keval,wkcoef,wkpows,N)-1;
        gn=eval_J_NDk_Nmodes(keval,Jcoef,Jpows,N)';
        
        
        
        if nargin>=7
            if toc(timer)>timeout
                error(['timeout',num2str(toc(timer))])
            end
        end       
    else
        n=[];
        gn=[];
    end
    
end


