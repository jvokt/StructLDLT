function As = Order4Symmetrize(A,n)
  A = reshape(A, n,n, n,n);
  A = A + permute(A, [2; 1; 4; 3]);
  As = reshape(A, n^2, n^2)/2;
end
