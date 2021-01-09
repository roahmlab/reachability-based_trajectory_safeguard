function J_out = eval_J_NDk_Nmodes(k,Jcoeffs,Jpowers,N)
% J_out = eval_J_2Dk_Nmodes(k,Jcoeffs,Jpowers,N)
%
% Given N modes and k nK-by-1, find the Jacobian
    
Ntotal = cumsum(N) ;
Nk = length(k);
Npowers = size(Jpowers,1);
kmat = repmat(k',Npowers,1);

if Ntotal(end)>0
    
    J_out = zeros(Ntotal(end),Nk) ;

    if N(1)>0
    
        for j=1:Nk
            J_out(1:N(1),j) = Jcoeffs{1}(N(1)*(j-1)+1:N(1)*j,:)*(prod(kmat.^Jpowers{1}(:,Nk*(j-1)+1:Nk*j),2)) ;
        end
    end
    
    for i = 2:length(N)
        if N(i)>0
            for j=1:Nk
                J_out((Ntotal(i-1)+1):Ntotal(i),j) = Jcoeffs{i}(N(i)*(j-1)+1:N(i)*j,:)*(prod(kmat.^Jpowers{i}(:,Nk*(j-1)+1:Nk*j),2)) ;
            end
         
        end
    end

else
   J_out=[]; 
end

end