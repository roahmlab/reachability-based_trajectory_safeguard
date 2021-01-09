function [Areduced,breduced] = reduce_row_rank(A,b)
Atranspose = A.';
[~,RB] = rref(Atranspose);
Areduced = Atranspose(:,RB).';
breduced=b(RB);
end

