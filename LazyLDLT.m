function [L,d,piv] = LazyLDLT(F, tol, n, A)
% Computes the nxn unit lower triangular L, n-vector d, 
% and n-vector piv where A(piv,piv) = L*diag(d)*L' using a rank-revealing,
% lazy-evaluation, symmetric-pivoting algorithm
% d(i) is zero if d(i) <= n*d(1)*tol
% F(i,j) is a function which returns the (i,j) entry of nxn matrix A

v = zeros(n,1);
w = zeros(n,1);
d = zeros(n,1);
L = eye(n,n);
piv = 1:n;

% Get the diagonal elements
for i=1:n
    d(i) = F(i,i,A);
end

j = 1;
% While the most recently computed value of d is greater than the threshold
while j <= n && (j == 1 || d(j-1) > n*d(1)*tol)
    % Search d(j:n) for largest diagonal element
    [alpha,idx] = max(d(j:n));
    k = idx+j-1;
    
    % Swap dk and dj
    d([k j]) = d([j k]);
    
    % Update piv
    piv([k j]) = piv([j k]);
    
    % Pivot L
    L([k j],:) = L([j k],:);
    L(:,[k j]) = L(:,[j k]);
    
    % Compute w, the next subcolumn of permuted A
    for i=j+1:n
        w(i) = F(piv(i),piv(j),A);
    end
    
    % Compute d(j) and L(:,j)
    for i=1:j-1
        v(i)= L(j,i)*d(i);
    end
    d(j) = alpha - L(j,1:j-1)*v(1:j-1);
    L(j+1:n,j) = (w(j+1:n) - L(j+1:n,1:j-1)*v(1:j-1))/d(j);
    
    j = j+1;
end

if (d(j-1) > n*d(1)*tol)
    j = j-1;
else
    j = j-2;
end
L = L(:,1:j);
d = d(1:j);