function ShowApproxQuantumSum()
n = 10;
[approxtime naivetime r error] = CompareQuantumSum(n)

function [approxtime naivetime r error] = CompareQuantumSum(n)
% Tests ApproxQuantumSum
% Software dependencies: rand3sym, ApproxQuantumSum, MatToTen
%n = 2;
fprintf('Starting n=%d\n',n);
% A is a random positive definite matrix with three symmetries: block 
% symmetry, symmetric blocks, and perfect shuffle permutation symmetry
A = RandOrder4ERI(n);
% v is an input to ApproxQuantumSum
v = randn(n,1);
% delta determines when to stop the Quantum Summation
delta = 1e-13;
% tol determines when to stop LazyLDL
tol = 1e-13;

tic;
mu = ApproxQuantumSum(A,v,delta,tol);
approxtime = toc;

T = MatToTen(A,[1 3],[2 4],[n n n n]);

tic;
nu = NaiveQuantumSum(T,v);
naivetime = toc;

error = abs(mu - nu)/abs(nu);
r = rank(A);
% fprintf('Relative error: %d\n',error);

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

function y = F(i,j,A)
global numfevals;
numfevals = numfevals + 1;
y = A(i,j);

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
while j <= n && d(j) > n*d(1)*tol
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

if (j ~= 1 && d(j-1) > n*d(1)*tol)
    j = j-1;
elseif (j ~= 1)
    j = j-2;
end
L = L(:,1:j);
d = d(1:j);


function QTA = Qskew(Askew,n,~)
% Returns Qskew*A without matrix multiplication
% Software dependencies: StructReshapeQskewTimesVector
[~,n2] = size(Askew);
QTA = zeros(n^2,n2);
% Computes the resulting matrix column by column
for k=1:n2
    QTA(:,k) = reshape(StructReshapeQskewTimesVector(Askew(:,k),n),n^2,1);
end


function S = StructReshapeQskewTimesVector(v,n)
% Returns Qskew*v computed without matrix multiplication reshaped as n x n
% S is skew-symmetric since Qskew is an orthogonal basis for reshaped
% skewmetric matricies

S = zeros(n);
% The lower triangular part of S is enumerated as 
% the elements of v
S(find(tril(ones(n),-1))) = v/sqrt(2);
S = S - S';


function QTA = QskewT(A,n,nskew)
% Returns Qskew'*A without matrix multiplication
% Software dependencies: QskewTransposedTimesVector
[~,n2] = size(A);
QTA = zeros(nskew,n2);
% Computes the resulting matrix column by column
for k=1:n2
    QTA(:,k) = QskewTransposedTimesVector(A(:,k),n,nskew);
end

function w = QskewTransposedTimesVector(v,n,~)
% Returns Qsym'*v without matrix multiplication
w = zeros(n*(n+1)/2,1);
Ak = reshape(v,n,n);
% The first n elements of w are the diagonal elements of Ak
%w(1:n) = diag(Ak);
AkT = Ak';
% v1 is a vector of the the lower triangular elements of Ak
v1 = Ak(find(tril(ones(n),-1)));
% v2 is a vector of the the upper triangular elements of Ak
v2 = AkT(find(tril(ones(n),-1)));
% w is a linear combo of v1 and v2
w = (v1 - v2)/sqrt(2);


function QTA = Qsym(Asym,n,~)
% Returns Qsym*Asym without matrix multiplication
% Software dependencies: StructReshapeQsymTimesVector
[~,n2] = size(Asym);
QTA = zeros(n^2,n2);
% Computes the resulting matrix column by column
for k=1:n2
    QTA(:,k) = reshape(StructReshapeQsymTimesVector(Asym(:,k),n),n^2,1);
end

function S = StructReshapeQsymTimesVector(v,n)
% Returns Qsym*v computed without matrix multiplication reshaped as n x n
% S is symmetric since Qsym is an orthogonal basis for reshaped
% symmetric matricies

% The diagonal elements of S are the first n elements of v
D = diag(v(1:n));
S = zeros(n);
% The lower triangular part of S is enumerated as 
% the remaining elements of v
S(find(tril(ones(n),-1))) = v(n+1:end)/sqrt(2);
S = S + D + S';



function QTA = QsymT(A,n,nsym)
% Returns Qsym'*A without matrix multiplication
% Software dependencies: QsymTransposedTimesVector
[~,n2] = size(A);
QTA = zeros(nsym,n2);
% Computes the resulting matrix column by column
for k=1:n2
    QTA(:,k) = QsymTransposedTimesVector(A(:,k),n,nsym);
end


function w = QsymTransposedTimesVector(v,n,nsym)
% Returns Qsym'*v without matrix multiplication

Ak = reshape(v,n,n);
% The first n elements of w are the diagonal elements of Ak
w(1:n) = diag(Ak);
AkT = Ak';
% v1 is a vector of the the lower triangular elements of Ak
v1 = Ak(find(tril(ones(n),-1)));
% v2 is a vector of the the upper triangular elements of Ak
v2 = AkT(find(tril(ones(n),-1)));
% The remaining elements of w are sqrt(2) times the average of v1 and v2
w(n+1:nsym) = (v1 + v2)/sqrt(2);


function [B, sig, C, tau] = StructLDLT(A,tol)
% A = sum_k sig(k) kron(Bk,Bk) + sum_k tau(k) kron(Ck,Ck)
% where Bk = reshape(B(:,k),n,n), Ck = reshape(C(:,k),n,n)
% Software dependencies: LazyLDLT, QsymT, StructReshapeQsymTimesVector,
% QskewT, StructReshapeQskewTimesVector

[n2,~] = size(A);
n = sqrt(n2);

nsym = n*(n+1)/2;
nskew = nsym - n;

if nargout < 3
    % Computing Asym = Qsym'*A*Qsym
    % QsymTA = Qsym'*A
    QsymTA = QsymT(A,n,nsym);
    % Asym = (Qsym'*(Qsym'*A)')'
    Asym = QsymT(QsymTA',n,nsym)';
%{
	Q = full(Q3sym(7));
	Qsym = Q(:,1:28);
	Asym = Qsym'*A*Qsym;
%}
    % Lazy-evaluation LDLT factorization
    [Lsym,sig,psym] = LazyLDLT(@F, tol, nsym, Asym);
%{
	Lsym = chol(Asym,'lower');
	AsymRank = rank(Asym)
	sig = diag(Lsym);
	Lsym = Lsym*diag(1./sqrt(sig));
	LDLTerror = norm(Lsym*diag(sig)*Lsym'-Asym)
%}	
    % rsym is the revealed rank of Asym
    rsym = length(sig);
%	psym = 1:rsym;
    v = zeros(rsym,1);
    B = cell(rsym,1);
%	QsymLsym = Qsym*Lsym;
    for k=1:rsym
%	B{k} = reshape(QsymLsym(:,k),n,n);
        % v = P'*Lsym(:,k)
        v(psym) = Lsym(:,k);
        % Reshapes Qsym*v computed without matrix multiplication
        B{k} = StructReshapeQsymTimesVector(v,n);
%	B{k} = reshape(Qsym*v,n,n);
    end
else
%	disp('Computing the entire factorization')
    % Computing Asym = Qsym'*A*Qsym
    % QsymTA = Qsym'*A
    QsymTA = QsymT(A,n,nsym);
    % Asym = (Qsym'*(Qsym'*A)')'
    Asym = QsymT(QsymTA',n,nsym)';
    % Lazy-evaluation LDLT factorization
%	disp('Lazy-evaluation LDLT factorization')
    [Lsym,sig,psym] = LazyLDLT(@F, tol, nsym, Asym);
    % rsym is the revealed rank of Asym
    rsym = length(sig);
    v = zeros(rsym,1);
    B = cell(rsym);
%	disp('computing B factors')
    for k=1:rsym
        % v = P'*Lsym(:,k)
        v(psym) = Lsym(:,k);
        % Reshapes Qsym*v computed without matrix multiplication
        B{k} = StructReshapeQsymTimesVector(v,n);
    end
    % Computing Askew = Qskew'*A*Qskew
    % QskewTA = Qskew'*A
    QskewTA = QskewT(A,n,nskew);
	normQskewTA = norm(QskewTA);
    % Askew = (Qskew'*(Qskew'*A)')'
    Askew = QskewT(QskewTA',n,nskew)';
	normAskew = norm(Askew);
    % Lazy-evaluation LDLT factorization	
%	disp('Lazy-evaluation LDLT factorization')
    [Lskew,tau,pskew] = LazyLDLT(@F, tol, nskew, Askew);
    % rskew is the revealed rank of Askew
    rskew = length(tau);
    v = zeros(rskew,1);
    C = cell(rskew);
%	disp('computing D factors')
    for k=1:rskew
        % v = P'*Lskew(:,k)
        v(pskew) = Lskew(:,k);
        % Reshapes Qskew*v computed without matrix multiplication
       C{k} = StructReshapeQskewTimesVector(v,n);
    end
end

function nu = NaiveQuantumSum(T,v)
nu = 0;
n = length(v);
for i=1:n
    for j=1:n
        for k=1:n
            for l=1:n
                nu = nu + T(i,j,k,l)*v(i)*v(j)*v(k)*v(l);
            end
        end
    end
end

function A = RandOrder4ERI(n)

A = randn(n,n^2);
A = A'*A;
A = Order4Symmetrize(A,n);

function As = Order4Symmetrize(A,n)
  A = reshape(A, n,n, n,n);
  A = A + permute(A, [2; 1; 4; 3]);
  As = reshape(A, n^2, n^2)/2;

function A = MatToTen(A_rxc,r,c,n)
% function A = MatToTen(A_rxc,r,c,n)
% A_rxc is a matrix.
% n is a positive integer d-vector
% r and c are integer vectors with the property that [r c] is a 
% permutation of 1:d.
% A_rxc is a matrix with prod(n(r)) rows and prod(n(c)) columns.
% A is an n(1)-by- ... -by-n(d) tensor where 
%       A_rxc(col(i(r),n(r)),col(i(c),n(c)) = A(i) 
% for all i that satisfy 1<=i<=n
% GV4: Section 12.4.6

% The key parameters..
d = length(n);
i = zeros(1,d);
nr = n(r);
nc = n(c);
[nRows,nCols] = size(A_rxc);
% Compute Vec(A)...
A = zeros(prod(n),1);
for p=1:nRows
    for q=1:nCols
        ir = inv_col(p,nr);
        ic = inv_col(q,nc);
        i(r) = ir;
        i(c) = ic;
        k = col(i,n);
        A(k) = A_rxc(p,q);
    end  
end
% Reshape into an order-d tensor...
A = reshape(A,n);

function A_rxc = TenToMat(A,r,c)
% function A_rxc = TenToMat(A,r,c)
% A is an n(1)x...xn(d) tensor.
% r and c are integer vectors with the property that [r c] is a 
% permutation of 1:d.
% A_rxc is a matrix with prod(n(r)) rows and prod(n(c)) columns defined by
% A_rxc(col(i(r),n(r)),col(i(c),n(c)) = A(i) for all i that satisfy
% 1<=i<=n. That is, A_rxc is the r-by-c unfolding of A.
% GV4: Section 12.4.6

% The key parameters...
n = size(A);
d = length(n);
i = zeros(1,d);
nr = n(r);
nRows = prod(nr);
nc = n(c);
nCols = prod(nc);
A_rxc = zeros(nRows,nCols);
a = Vec(A);
% Set up the unfolding...
for p=1:nRows
    for q=1:nCols
        i(r) = inv_col(p,n(r));
        i(c) = inv_col(q,n(c));
        A_rxc(p,q) = a(col(i,n));
    end
end

function k = col(i,n)
% function k = col(i,n)
% i and n are integer d-vectors with 1<=i<=n
% k = i(1) + (i(2)-1)n(1) + (i(3)-1)n(1)n(2) + ... + (i(d)-1)n(1)...n(d-1)
% GVL4: Section 12.4.5
d = length(n);
k = i(d)-1;
for j=d-1:-1:1
   k = n(j)*k + i(j)-1;
end
k = k+1;

function i = inv_col(k,n)
% function i = inv_col(k,n)
% n is a d-vector with positive integer entries
% k is an integer that satisfies 1<=k<=prod(n)
% i is a d-vector with positive integer entries so
% k = i(1) + (i(2)-1)n(1) + (i(3)-1)n(1)n(2) + ... + (i(d)-1)n(1)...n(d-1)
% GVL4: Section 12.4.5
d = length(n);
i = zeros(d,1);
N = prod(n);
k = k-1;
for j=d:-1:1
    N = N/n(j);
    f = floor(k/N);
    k = k-N*f;
    i(j) = f+1;   
end

function Q = Q3sym(n)
% Q is the n^2xn^2 orthonormal matrix defined in the UCLA pdf

Q = sparse(n^2,n^2);
for i=1:n
    Q(1+(n+1)*(i-1),i) = 1;
end

k = n+1;
p = n*(n+1)/2+1;
for i=1:n
  for j=i+1:n
    Q(i+(j-1)*n,k) = 1/sqrt(2);
    Q(j+(i-1)*n,k) = 1/sqrt(2);
    Q(i+(j-1)*n,p) = -1/sqrt(2);
    Q(j+(i-1)*n,p) = 1/sqrt(2);
    k = k+1;
    p = p+1;
  end
end

function v = PerfShuff(p,r)
% function v = PerfShuff(p,r)
% Vector representation of the perfect shuffle permutation.
% GVL4: Setion 1.2.11
% p and r are positive integers. If n = pr then v is a length-n row
%   vector with the property that if x is a column n-vector and
%   y = x(v), then y = I(v,:)*x = P_{p,r}* x is the mod-p perfect shuffle
%   of x.
n = p*r;
v = [];
for k=1:r
    v = [v k:r:n];
end