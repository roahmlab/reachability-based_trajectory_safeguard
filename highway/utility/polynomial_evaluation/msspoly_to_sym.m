function q= msspoly_to_sym(p,x,xsym)
% this function takes in an msspoly and returns it in matlab's sym data
% type

%inputs
%p: n by k msspoly
%x: nvars by 1 msspoly containing variables in p
%xsym: nvars by 1 symbolic variables matching x

%outputs
%q: n by k symbolic variable s.t. q(xsym)==p(x)

if all(size(x)==size(xsym))
    [vars,pow,M]=decomp(p);
    
    var_idxs = NaN(size(vars));
    for i=1:length(vars)
        var_idxs(i) = find(arrayfun(@(xin)isequal(vars(i),xin),x));
    end
    
    q=sym(zeros(size(p)));
    
    for i=1:size(M,1)
        
        for k=1:size(M,2)
            term = sym(1);
            for j=1:size(pow,2)
                term = term * xsym(var_idxs(j))^pow(k,j);
            end
            q(i)=q(i)+M(i,k)*term;
        end
    end
else
    error('size of x and xsym must match!')
end

end