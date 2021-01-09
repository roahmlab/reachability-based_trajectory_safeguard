function W_out = eval_w_many_points(K,wc,wp)
% W_out = eval_w_many_points(K,wc,wp)
%
% Evaluate a polynomial w(k), expressed as a coefficient matrix (wc) and a
% powers matrix (wp), at the points K (2-by-N). This is the same as the
% function eval_w, but can handle a 2-by-N input. However, calling this
% function on a single point is 3 to 4 times as slow as eval_w.
%
% The output W_out is N-by-1.

    k1k2eval = repmat(K(:)',size(wp,1),1).^repmat(wp,1,size(K,2)) ;
    kprod = k1k2eval(:,1:2:end).*k1k2eval(:,2:2:end) ;
    W_out = wc*kprod ;
end