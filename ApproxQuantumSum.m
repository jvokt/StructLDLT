function [mu,nu] = ApproxQuantumSum(A,v,delta,tol)
% A is a uniformly blocked n^2xn^2 matrix with 3-symmetry: 
% nxn symmetric blocks, block symmetry, and perfect shuffle symmetry
% Computes mu = sum_i sum_j sum_k sum_l A(i,j,k,l) v(i) v(j) v(k) v(l)
% where v is a column-n vector and 
% A is a reshaped fourth-order tensor with the following symmetries: 
% A(i,j,k,l) = A(j,i,k,l) = A(i,j,l,k) = A(k,l,i,j)
% A is assumed to be unfolded as A(i,j,k,l) = [A_k,l]_i,j
% Tensor is positive semidefinite, relatively low-rank
% Software dependencies: StructLDLT

if nargout < 2
    % Computes the Kronecker Product SVD using lazy-evaluation LDL^T
    [B, sig] = StructLDLT(A,tol);
    mu = 0;
    s = delta + 1;
    % Loops while mu is changing by more than delta each iteration
    k = 1;
    while s > delta && k <= length(sig)
        s = sig(k)*(v'*(B{k}*v))^2;
        mu = mu + s;
        k = k + 1;
    end
else
    % Computes the Kronecker Product SVD using lazy-evaluation LDL^T
    [B, sig, C, tau] = StructLDLT(A,tol);    
    mu = 0;
    s = delta + 1;
    % Loops while mu is changing by more than delta each iteration
    k = 1;
    while s > delta && k <= length(sig)
        s = sig(k)*(v'*(B{k}*v))^2;
        mu = mu + s;
        k = k + 1;
    end
    nu = 0;
    s = delta + 1;
    % Loops while nu is changing by more than delta each iteration
    k = 1;
    while s > delta && k <= length(tau)
        s = tau(k)*(v'*(C{k}*v))^2;
        nu = nu + s;
        k = k + 1;
    end
end