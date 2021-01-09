function samples = rand_semialg(h,box,N,x)
    
%create random numbers in R^n that lie in msspoly semialegraic sets
%input: h n-by-1 msspoly of polynomials defining sets (h{i}(x)>0) or cell
%       array of function handles
%       box n-by-2 box to generate random numbers in (this should be larger
%       than 0 super-level set of interest in semi-algebraic set
%       N number of samples to generate
%       x n-by-1 msspoly variables fo h (if msspoly)

N_h = length(h);
n = size(box,1);


% convert msspoly to function handles
if iscell(h)
    h_fun = h;
else
    h_fun = cell(N_h,1);
    for i = 1:length(h)
        h_fun{i} = msspoly_to_fun(h(i),{x});
    end
end


samples = randRange(box(:,1),box(:,2),[],[],1,N);
for i = 1:length(h)
    L_fail = h_fun{i}(samples) < 0;
    samples(:,L_fail) = NaN(n,sum(L_fail));
end

while any(isnan(samples(1,:)))
    L_nan = isnan(samples(1,:));
    samples(:,L_nan) = randRange(box(:,1),box(:,2),[],[],1,sum(L_nan));
    for i = 1:length(h)
        L_fail = h_fun{i}(samples) < 0;
        samples(:,L_fail) = NaN(n,sum(L_fail));
    end
end

end