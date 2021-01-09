function J_out = eval_J_2Dk(k,Jcoeffs,Jpowers,N)
    J_out = zeros(N,2) ;
    Npowers = size(Jpowers,1) ;
    kmat = repmat(k',Npowers,1) ;
    J_out(:,1) = Jcoeffs(1:N,:)*(prod(kmat.^Jpowers(:,1:2),2)) ;
    J_out(:,2) = Jcoeffs((N+1):(2*N),:)*(prod(kmat.^Jpowers(:,3:4),2)) ;
end