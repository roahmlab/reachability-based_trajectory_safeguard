function J_out = eval_J(k,Jcoeffs,Jpowers,N)
% J_out = eval_J(k,Jcoeffs,Jpowers,N)
%
% Given a Jacobian matrix represented as a coefficients matrix and a powers
% matrix (see diff_wk_wrt_k), evaluate the jacobian at k (Nk-by-1)

% NEW VERSION, NO TRANSPOSE OF JCOEFFS
    Nk=length(k);
    
    J_out = zeros(N,Nk) ;
% 
    for i=1:Nk
        J_out(:,i) = Jcoeffs(1+N*(i-1):N*i,:)*(prod(repmat(k',size(Jpowers,1),1).^Jpowers(:,(i-1)*Nk+1:Nk*i),2)) ;
    end

% SEAN'S OLD VERSION
%     Nk=length(k);
%     
%     J_out = zeros(N,Nk) ;
% % 
%     for i=1:Nk
%         J_out(:,i) = Jcoeffs(:,1+N*(i-1):N*i)'*(prod(repmat(k',size(Jpowers,1),1).^Jpowers(:,(i-1)*Nk+1:Nk*i),2)) ;
%     end

% SHREYAS' OLD VERSION
%     J_out(:,1) = Jcoeffs(:,1:N)'*(prod(repmat(k',size(Jpowers,1),1).^Jpowers(:,1:2),2)) ;
%     J_out(:,2) = Jcoeffs(:,N+1:end)'*(prod(repmat(k',size(Jpowers,1),1).^Jpowers(:,3:4),2)) ;
end