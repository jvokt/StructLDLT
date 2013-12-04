function A = Rand3SymAlt(n,d1,d2)
% A is a random n^2xn^2 indefinite 3-symmetric matrix
% with Kronecker lambda(A) = lambda(diag(d1)) union lambda(diag(d2))

nsym = n*(n+1)/2;
nskew = nsym - n;

if nargin < 2
    d1 = rand(nsym,1);
    d2 = rand(nskew,1);
end

A = zeros(n^2,n^2);

for i=1:length(d1)
    B = randn(n,n);
    B = B + B';
    A = A + d1(i)*(kron(B,B));
end

for i=1:length(d2)
    C = randn(n,n);
    C = C + C';
    A = A + d2(i)*(kron(C,C));
end