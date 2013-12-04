function e = Order6symmetric(A)
% Permutations on three indices
G = [1, 1, 2, 2, 3, 3;
     2, 3, 1, 3, 1, 2;
     3, 2, 3, 1, 2, 1];

sym_err = zeros(6,1);

for i=1:6
    idx = [G(:,i),G(:,i)+3];
    T = reshape(A,[n n n n n n]);
    sym_err(i) = norm(A-reshape(permute(T,idx),n^3,n^3));
end

e = norm(sym_err);