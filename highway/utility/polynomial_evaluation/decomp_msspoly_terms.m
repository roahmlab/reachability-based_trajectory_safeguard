    function [monomials,M]=decomp_msspoly_terms(p,terms)
         [varlist,~,~]=decomp(p);
         
         [~,varlistid]=isfree(varlist);
         [~,kid]=isfree(terms);
         
         L=~ismember(varlistid,kid);
        
         [x,p,M]=decomp(p,varlist(L));
         
         monomials=msspoly(zeros(length(M),1));
         for j=1:length(M)
             try
                monomials(j)=prod(x'.^p(j,:));
             catch
                 xtemp=x'.^p(j,:);
                 monomials(j)=1;
                 for k=1:length(xtemp)
                     monomials(j)=monomials(j)*xtemp(k);
                 end
             end
         end
        M=M';
    end