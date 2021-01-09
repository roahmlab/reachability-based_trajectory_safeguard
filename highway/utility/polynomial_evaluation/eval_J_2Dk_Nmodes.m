function J_out = eval_J_2Dk_Nmodes(k,Jcoeffs,Jpowers,N)
% J_out = eval_J_2Dk_Nmodes(k,Jcoeffs,Jpowers,N)
%
% Given N modes and k 2-by-1, find the Jacobian
    
Ntotal = cumsum(N) ;
if Ntotal(end)>0
    
    J_out = zeros(Ntotal(end),2) ;

    if N(1)>0

    J_out(1:N(1),1) = Jcoeffs{1}(1:N(1),:)*(prod(k(1).^Jpowers{1}(:,1:2),2)) ;
    J_out(1:N(1),2) = Jcoeffs{1}(N(1)+1:end,:)*(prod(k(2).^Jpowers{1}(:,3:4),2)) ;
    end
    
    for i = 2:length(N)
        if N(i)>0
            J_out((Ntotal(i-1)+1):Ntotal(i),1) = Jcoeffs{i}(1:N(i),:)*(prod(k(1).^Jpowers{i}(:,1:2),2)) ;
            J_out((Ntotal(i-1)+1):Ntotal(i),2) = Jcoeffs{i}((N(i)+1):end,:)*(prod(k(2).^Jpowers{i}(:,3:4),2)) ;
        end
    end

else
   J_out=[]; 
end

end