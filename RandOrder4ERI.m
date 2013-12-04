function A = RandOrder4ERI(n)

A = randn(n^2,n^2);
A = A'*A;
A = Order4Symmetrize(A,n);
