function As = Order6symmetrize(A,n)
% Permutations on three indices
G = [1, 1, 2, 2, 3, 3;
     2, 3, 1, 3, 1, 2;
     3, 2, 3, 1, 2, 1];

  A = reshape(A, n,n,n, n,n,n);
  As = A;
  idx1 = [1;2;3];
  idx2 = [4;5;6];
  for k = 2:6
    idx = [idx1(G(:,k)); idx2(G(:,k))];
    As = As + permute(A, idx);
  end
  As = reshape(As, n^3, n^3)/6;
end