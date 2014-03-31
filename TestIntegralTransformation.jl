using PyPlot

function TestIntegralTransformation()
# Calls CompareIntegralTransformation for different values of n

maxn = 25;
comparedtimes = zeros(maxn)
ourtimes = zeros(maxn)
errors = zeros(maxn)
ranks = zeros(maxn)
fevals = zeros(maxn)
global numfevals
for n=2:maxn
    c, o, e, r = CompareIntegralTransformation(n)
    comparedtimes[n] = c
    ourtimes[n] = o
    errors[n] = e
    ranks[n] = r
    fevals[n] = numfevals
end
figure()
plot(1:maxn,ourtimes,color="red",label="StructLDLT")
plot(1:maxn,comparedtimes,color="blue",label="Full SVD")
legend()
xlabel("n=2 to 25")
ylabel("Time in seconds")
title("Time in seconds to tranform all integrals")
figure()
plot(1:maxn,-log10(errors));
xlabel("n=2 to 25")
ylabel("-Log(error)")
title("Frobenius norm error in transformed matrix")
figure()
plot(1:maxn,ranks);
xlabel("n=2 to 25")
ylabel("Numerical rank")
title("Rank of the ERI matrix")
figure()
plot(1:maxn,fevals);
xlabel("n=2 to 25")
ylabel("F evals")
title("Number of two-electron integral evaluations")
    
end

function CompareIntegralTransformation(n)
# Compares SVDIntegralTransformation with FastIntegralTransformation

# Random n^2 by n^2 ERI matrix of rank 2n with correct symmetries
# Unfolded as A_ij,kl=(ij|kl)
A = RandERIMatrix(n);
r = rank(A)
# Coefficient matrix from Hartree-Fock method
C = randn(n,n);
# tol determines when to stop LazyLDLT
tol = 1e-13;
# Reshapes ERI tensor as T[i,j,k,l]=(ij|kl)
T = reshape(A,n,n,n,n)

tic()
M2,D,Z = SVDIntegralTransformation(A,C,tol);
naivetime = toc();

tic()
M1,sig,Z = FastIntegralTransformation(A,C,tol);
fasttime = toc();

error = normfro(M1-M2);
return naivetime, fasttime, error, r
end

function RandERIMatrix(n)
# Returns a random positive definite ERI matrix with the appropriate symmetries
# Rank is hard coded as 2n
    r = 2n
    A = zeros(n^2,n^2)
    for i=1:r
        B = randn(n,n)
        B = (B + B')/2
        A += vec(B)*vec(B)'
    end
    return A
end

function SortOfNaiveIntegralTransformation(A,C)
return kron(C,C)'*A*kron(C,C)
end

function ReallyNaiveIntegralTransformation(T,C)
n,n,n,n = size(T)
M = zeros(n,n,n,n)
for p=1:n
    for q=1:n
        for r=1:n
            for s=1:n
                for mu=1:n
                    for nu=1:n
                        for lam=1:n
                            for sig=1:n
                                M[p,q,r,s] += C[p,mu]*C[q,nu]*C[r,lam]*C[s,sig]*T[mu,nu,lam,sig]
                            end
                        end
                    end
                end
            end
        end
    end
end
return reshape(M,n^2,n^2)
end

function NaiveIntegralTransformation(T,C)
    n,n,n,n = size(T)

    MO = zeros(n,n,n,n)
    temp = zeros(n,n,n,n)
    temp2 = zeros(n,n,n,n)
    temp3 = zeros(n,n,n,n)
    for mu = 1:n
	for i = 1:n
	    temp[mu,:,:,:] += C[i,mu]*T[i,:,:,:]
	for nu = 1:n
	    for j = 1:n
		temp2[mu,nu,:,:] += C[j,nu]*temp[mu,j,:,:]
	    for lam = 1:n
		for k = 1:n
		    temp3[mu,nu,lam,:] += C[k,lam]*temp2[mu,nu,k,:]
		for sig = 1:n
		    for l = 1:n
			MO[mu,nu,lam,sig] += C[l,sig]*temp3[mu,nu,lam,l]
                    end
                end
                end
            end
            end
        end
        end
    end
    return reshape(MO,n^2,n^2)
end

function SVDIntegralTransformation(A,C,tol)
# Gets the SVD of A and returns 
# 1. The transformed matrix M
# 2. The eigenvalues of A
# 3. A cell array of the Kronecker Product SVD factors
    n2,n2 = size(A)
    n = convert(Int64,sqrt(n2))
    (V,D) = svd(A);
    Z = cell(n2)
    M = zeros(n2,n2)
    Mt = zeros(n2,n2)
    for k=1:n2
        v = V[:,k]
        B = reshape(v,n,n)        
        Z[k] = (C'*(B*C))
        v = vec(Z[k])
        M = M + D[k]*v*v'
    end
return M,D,Z
end

function FastIntegralTransformation(A,C,tol)
# A is an n^2 by n^2 ERI matrix unfolded as A_ij,kl=(ij|kl)
# C is the coefficent matrix from the Hartree Fock method
# Transforms A in Atomic Orbital basis to M in Molecular Orbital basis
    n2,n2 = size(A)
    # Computes the Kronecker Product SVD using lazy-evaluation LDL^T
    (B, sig) = StructLDLT(A,tol);
    rsym = length(sig)
    Z = cell(rsym)
    M = zeros(n2,n2)
    Mt = zeros(n2,n2)
    for k=1:rsym
        Z[k] = (C'*(B[k]*C))
        v = vec(Z[k])
        M = M + sig[k]*v*v'
    end
return M,sig,Z
end

function StructLDLT(A,tol)
# A = sum_k sig(k) kron(B[k],B[k])

(n2,x) = size(A);
n = convert(Int64,sqrt(n2));

nsym = convert(Int64,n*(n+1)/2);
nskew = nsym - n;

    # Computing Asym = Qsym'*A*Qsym
    # QsymTA = Qsym'*A
    QsymTA = QsymT(A,n,nsym);
    # Asym = (Qsym'*(Qsym'*A)')'
    Asym = QsymT(QsymTA',n,nsym)';
    # Lazy-evaluation LDLT factorization
    (Lsym,sig,psym) = LazyLDLT(F, tol, nsym, Asym);
    # rsym is the revealed rank of Asym
    rsym = length(sig);
    v = zeros(nsym);
    B = cell(rsym);
    for k=1:rsym
        # v = P'*Lsym(:,k)
        v[psym] = Lsym[:,k];
        # Reshapes Qsym*v computed without matrix multiplication
        B[k] = StructReshapeQsymTimesVector(v,n);
    end
return (B, sig)
end

function QsymT(A,n,nsym)
# Returns Qsym'*A without matrix multiplication
(x,n2) = size(A);
QTA = zeros(nsym,n2)
# Computes the resulting matrix column by column
for k=1:n2
    QTA[:,k] = QsymTransposedTimesVector(A[:,k],n,nsym);
end
return QTA
end

function QsymTransposedTimesVector(v,n,nsym)
# Returns Qsym'*v without matrix multiplication
w = zeros(nsym);
Ak = reshape(v,n,n);
nskew = n^2-nsym;
v1 = zeros(nskew);
v2 = zeros(nskew);
# The first n elements of w are the diagonal elements of Ak
d = diag(Ak);
w[1:n] = d;
AkT = Ak';
# v1 is a vector of the the lower triangular elements of Ak
for j=1:n-1
    v1[n*(j-1)-j*(j-1)/2+1:n*j-j*(j+1)/2] = Ak[j+1:n,j];
end
# v2 is a vector of the the upper triangular elements of Ak
for j=1:n-1
    v2[n*(j-1)-j*(j-1)/2+1:n*j-j*(j+1)/2] = AkT[j+1:n,j];
end
# The remaining elements of w are sqrt(2) times the average of v1 and v2
w[n+1:nsym] = (v1 + v2)/sqrt(2);
return w
end

function F(i,j,A)
global numfevals;
numfevals = numfevals + 1;
return A[i,j];
end

function LazyLDLT(F, tol, n, A)
# Computes the nxn unit lower triangular L, n-vector d, 
# and n-vector piv where A(piv,piv) = L*diag(d)*L' using a rank-revealing,
# lazy-evaluation, symmetric-pivoting algorithm
# d(i) is zero if d(i) <= n*d(1)*tol
# F(i,j) is a function which returns the (i,j) entry of nxn matrix A
global numfevals = 0;
v = zeros(n);
w = zeros(n);
d = zeros(n);
L = eye(n,n);
piv = [i for i in 1:n];
thresh = n*d[1]*tol
# Get the diagonal elements
for i=1:n
    d[i] = F(i,i,A);
end
j = 1;
crit = d[1]
while j <= n && crit > thresh
    # Search d(j:n) for largest diagonal element
    if j == n
        alpha = d[j]
        idx = 1
    else
        dvec = d[j:n]
        (alpha,idx) = maximum(dvec), indmax(dvec)
    end
    k = idx+j-1;
    # Swap dk and dj
    # d([k j]) = d([j k]);
    d[k], d[j] = d[j], d[k];
    # Update piv
    # piv([k j]) = piv([j k]);
    piv[k], piv[j] = piv[j], piv[k];
    # Pivot L
    # L([k j],:) = L([j k],:);
    L[k,:], L[j,:] = L[j,:], L[k,:];
    # L(:,[k j]) = L(:,[j k]);
    L[:,k], L[:,j] = L[:,j], L[:,k];
    # Compute w, the next subcolumn of permuted A
    for i=j+1:n
        w[i] = F(piv[i],piv[j],A);
    end
    # Compute d(j) and L(:,j)
    for i=1:j-1
        v[i]= L[j,i]*d[i];
    end
    c = L[j,1:j-1]*v[1:j-1]
    d[j] = alpha - c[1]
    L[j+1:n,j] = (w[j+1:n] - L[j+1:n,1:j-1]*v[1:j-1])/d[j];
    crit = d[j]
    j = j+1;
end
if (j != 1 && crit > thresh)
    j = j-1;
elseif (j != 1)
    j = j-2;
end
L = L[:,1:j];
d = d[1:j];
return (L,d,piv)
end

function Qsym(Asym,n,x)
# Returns Qsym*Asym without matrix multiplication
(x,n2) = size(Asym);
QTA = zeros(n^2,n2);
# Computes the resulting matrix column by column
for k=1:n2
    QTA[:,k] = reshape(StructReshapeQsymTimesVector(Asym[:,k],n),n^2,1);
end
return QTA
end

function StructReshapeQsymTimesVector(v,n)
# Returns Qsym*v computed without matrix multiplication reshaped as n x n
# S is symmetric since Qsym is an orthogonal basis for reshaped
# symmetric matricies

# The diagonal elements of S are the first n elements of v
D = diagm(v[1:n]);
S = zeros(n,n);
# The lower triangular part of S is enumerated as 
# the remaining elements of v
# S(find(tril(ones(n,n),-1))) = v[n+1:length(v)]/sqrt(2);
for j=1:n-1
    S[j+1:n,j] = v[n+n*(j-1)-j*(j-1)/2+1:n+n*j-j*(j+1)/2]/sqrt(2);
end
S = S + D + S';
return S
end