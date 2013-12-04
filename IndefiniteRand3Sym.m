function A = IndefiniteRand3Sym(n)
% Returns an random indefinite three symmetric matrix

A = zeros(n^2,n^2);
for i = 1:n
    for j = i:n
        B = randn(n,n);
        B = B + B';
        A((i-1)*n + 1:i*n, (j-1)*n + 1:j*n) = B;
        A((j-1)*n + 1:j*n, (i-1)*n + 1:i*n) = B;
    end
end
A = A + A(PerfShuff(n,n),PerfShuff(n,n));