function TestApproxQuantumSum()
% Tests ApproxQuantumSum
% Software dependencies: rand3sym, ApproxQuantumSum, MatToTen
n = 5;
% A is a random positive definite matrix with three symmetries: block 
% symmetry, symmetric blocks, and perfect shuffle permutation symmetry
A = rand3sym(n);
% v is an input to ApproxQuantumSum
v = randn(n,1);
% delta determines when to stop the Quantum Summation
delta = 1e-13;
% tol determines when to stop LazyLDL
tol = 1e-13;

mu = ApproxQuantumSum(A,v,delta,tol);

T = MatToTen(A,[1 3],[2 4],[n n n n]);

nu = 0;

for i=1:n
    for j=1:n
        for k=1:n
            for l=1:n
                nu = nu + T(i,j,k,l)*v(i)*v(j)*v(k)*v(l);
            end
        end
    end
end

error = abs(mu - nu)