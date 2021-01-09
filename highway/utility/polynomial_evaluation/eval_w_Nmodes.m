function n_out = eval_w_Nmodes(k,wcoeffs,wpowers,N)
% J_out = eval_J_2Dk_3modes(k,Jcoeffs,Jpowers,N)
%
% Given N modes and k 2-by-1, find the Jacobian


Ntotal = cumsum(N) ;
if Ntotal(end)>0
    
    n_out = zeros(Ntotal(end),1) ;
    
    if N(1)>0
        n_out(1:Ntotal(1),1) =  wcoeffs{1}*prod(repmat(k',size(wpowers{1},1),1).^wpowers{1},2) ;
    end
    
    for i=2:length(N)
        if N(i)>0
        n_out((Ntotal(i-1)+1):Ntotal(i),1) =  wcoeffs{i}*prod(repmat(k',size(wpowers{i},1),1).^wpowers{i},2) ;
        end
    end
  
else
   n_out=[]; 
end

end