function gout=remove_near_zero_terms(gin,tol,add)
if nargin<2
    tol=1e-6;
end
if nargin<3
    add='=';
end

tol=abs(tol);

gout=msspoly(zeros(numel(gin),1));
for j=1:length(gin)
[x,p,M]=decomp(gin(j));

L=abs(M)<tol;

Mnew=M;
Mnew(L)=0;

xnew=repmat(x',[size(p,1),1]);

varmat=xnew.^p;

gout(j)=0;

for i=1:size(p,1)
    try
    gout(j)=gout(j)+Mnew(i)*prod(varmat(i,:));
    catch
        pv=1;
        for k=1:length(varmat(i,:))
            pv=pv*varmat(i,k);
        end
         gout(j)=gout(j)+Mnew(i)*pv;
    end
end


if any(abs(Mnew)>0) || double(sum(abs(M(L))))>tol
if add=='+'
    gout(j)=gout(j)+tol*sum(L);
elseif add=='-'
    gout(j)=gout(j)-tol*sum(L);
end
else
warning('remove_near_zero terms: allterms below tolerance')
end

end
gout=reshape(gout,size(gin));
end