function J_out = eval_J_2Dk_4modes(k,Jcoeffs,Jpowers,N)
% J_out = eval_J_2Dk_3modes(k,Jcoeffs,Jpowers,N)
%
% Given 4 modes and k 2-by-1, find the Jacobian

% shreyas' new version
%     Nk = length(k);
    
    Ntotal = cumsum(N) ;
    if Ntotal(end)>0
    J_out = zeros(Ntotal(end),2) ;

%     k1vec = k(1).*ones(size(Jpowers,1),1) ;
%     k2vec = k(2).*ones(size(Jpowers,1),1) ;
    if N(1)>0

    J_out(1:N(1),1) = Jcoeffs{1}(1:N(1),:)*(prod(k(1).^Jpowers{1}(:,1:2),2)) ;
    J_out(1:N(1),2) = Jcoeffs{1}(N(1)+1:end,:)*(prod(k(2).^Jpowers{1}(:,3:4),2)) ;
    end
    
    if N(2)>0
    
    J_out((Ntotal(1)+1):Ntotal(2),1) = Jcoeffs{2}(1:N(2),:)*(prod(k(1).^Jpowers{2}(:,1:2),2)) ;
    J_out((Ntotal(1)+1):Ntotal(2),2) = Jcoeffs{2}((N(2)+1):end,:)*(prod(k(2).^Jpowers{2}(:,3:4),2)) ;
    end
    
    if N(3)>0
    J_out((Ntotal(2)+1):Ntotal(3),1) = Jcoeffs{3}(1:N(3),:)*(prod(k(1).^Jpowers{3}(:,1:2),2)) ;
    J_out((Ntotal(2)+1):Ntotal(3),2) = Jcoeffs{3}((N(3)+1):end,:)*(prod(k(2).^Jpowers{3}(:,3:4),2)) ;
    end
    
    if N(4)>0
    J_out((Ntotal(3)+1):Ntotal(4),1) = Jcoeffs{4}(1:N(4),:)*(prod(k(1).^Jpowers{4}(:,1:2),2)) ;
    J_out((Ntotal(3)+1):Ntotal(4),2) = Jcoeffs{4}((N(4)+1):end,:)*(prod(k(2).^Jpowers{4}(:,3:4),2)) ;
    end
    
    else
       J_out=[]; 
    end

end