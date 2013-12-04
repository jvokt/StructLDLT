function A = rand3sym(n,d1,d2)
% A is a random n^2xn^2 positive definite 3-symmetric matrix
% with lambda(A) = lambda(diag(d1)) union lambda(diag(d2))

nsym = n*(n+1)/2;
nskew = nsym - n;

if nargin < 2
    d1 = rand(nsym,1);
    d2 = rand(nskew,1);
end

[Q1,~] = qr(randn(nsym,nsym));
[Q2,~] = qr(randn(nskew,nskew));
Asym = Q1*diag(d1)*Q1';
Askew = Q2*diag(d2)*Q2';
QsymA = Qsym(Asym,n,nsym);
A1 = Qsym(QsymA',n,nsym)';
QskewA = Qskew(Askew,n,nskew);
A2 = Qskew(QskewA',n,nskew)';

A = A1 + A2;

for k=1:n
    for l=k:n
        I = n*(k-1)+1:n*k;
        J = n*(l-1)+1:n*l;
        A(I,J) = (A(I,J) + A(I,J)')/2;
        A(J,I) = (A(J,I) + A(J,I)')/2;
        A(I,J) = (A(I,J) + A(J,I))/2;
        A(J,I) = A(I,J);
    end
end
%A = A + A(PerfShuff(n,n),PerfShuff(n,n));