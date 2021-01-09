function n_out = eval_w_4modes(k,wcoeffs,wpowers,N)
% J_out = eval_J_2Dk_3modes(k,Jcoeffs,Jpowers,N)
%
% Given 3 modes and k 2-by-1, find the Jacobian

% shreyas' new version
%     Nk = length(k);
    
    Ntotal = cumsum(N) ;
    if Ntotal(end)>0
    n_out = zeros(Ntotal(end),1) ;

%     k1vec = k(1).*ones(size(Jpowers,1),1) ;
%     k2vec = k(2).*ones(size(Jpowers,1),1) ;
    if N(1)>0

    n_out(1:N(1),1) =  wcoeffs{1}*prod(repmat(k',size(wpowers{1},1),1).^wpowers{1},2) ;
    end
    
    if N(2)>0
    
    n_out((Ntotal(1)+1):Ntotal(2),1) =  wcoeffs{2}*prod(repmat(k',size(wpowers{2},1),1).^wpowers{2},2) ;
    end
    
    if N(3)>0
    n_out((Ntotal(2)+1):Ntotal(3)) =  wcoeffs{3}*prod(repmat(k',size(wpowers{3},1),1).^wpowers{3},2) ;
    end
    
    if N(4)>0
    n_out((Ntotal(3)+1):Ntotal(4)) =  wcoeffs{4}*prod(repmat(k',size(wpowers{4},1),1).^wpowers{4},2) ;
    end
    
    else
       n_out=[]; 
    end
% % sean's version
%     Nk=length(k);
%     
%     J_out = zeros(N,Nk) ;
% % 
%     for i=1:Nk
%     J_out(:,i) = Jcoeffs(:,1+N*(i-1):N*i)'*(prod(repmat(k',size(Jpowers,1),1).^Jpowers(:,(i-1)*Nk+1:Nk*i),2)) ;
%     end
end