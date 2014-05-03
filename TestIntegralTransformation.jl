using PyPlot

function TestIntegralTransformation()

    molecules = ["HF","NH3","H2O2","N2H4","C2H5OH"]
    num_basis_fns = [34, 48, 68, 82, 123]
    ranks = [200, 300, 400, 500, 750]
    nmols = 1# length(num_basis_fns)
    tol = 1e-6
    
    full_fact_time = zeros(nmols)
    full_trans_time = zeros(nmols)
    struct_fact_time = zeros(nmols)
    struct_trans_time = zeros(nmols)
    
    full_fact_space = zeros(nmols)
    full_trans_space = zeros(nmols)
    struct_fact_space = zeros(nmols)
    struct_trans_space = zeros(nmols)
    
    full_fact_error = zeros(nmols)
    full_trans_error = zeros(nmols)
    struct_fact_error = zeros(nmols)
    struct_trans_error = zeros(nmols)
    
    full_fact_fevals = zeros(nmols)
    struct_fact_fevals = zeros(nmols)
    
    for i=1:nmols
        molecule = molecules[i]
        n = num_basis_fns[i]
        r = ranks[i]
        A = RandERIMatrix(n,r)
        C = randn(n,n)
        
        # Full factorization
        tic()
        L, piv, full_fact_space[i], full_fact_fevals[i] = full_fact(A,tol)
        full_fact_time[i] = toc()
        full_fact_error[i] = norm(A[piv,piv] - L*L',1)
        
        # Full transformation
        tic()
        B, full_trans_space[i] = full_trans(L,piv,C)
        full_trans_time[i] = toc()
        full_trans_error[i] = norm(kron(C,C)*A*kron(C,C)' - B*B',1)
        
        # Structured factorization
        tic()
        L, piv, struct_fact_space[i], struct_fact_fevals[i] = struct_fact(A,tol)
        struct_fact_time[i] = toc()
        J = QsymTransform(L,piv,n)
        struct_fact_error[i] = norm(A - J*J',1)
        println("Structured factorization error: ",struct_fact_error[i])
        
        # Structured transformation
        tic()
        B, struct_trans_space[i] = struct_trans(L,piv,C)
        struct_trans_time[i] = toc()
        struct_trans_error[i] = norm(kron(C,C)*A*kron(C,C)' - B*B',1)
        println("Structured transformation error: ",struct_trans_error[i])
    end
    
    if true
        plotComparison(num_basis_fns[1:nmols],full_fact_time,"Full factorization",struct_fact_time,"Structured factorization","Number of basis functions","Time (seconds)","Time to factorize")
        plotComparison(num_basis_fns[1:nmols],full_fact_space,"Full factorization",struct_fact_space,"Structured factorization","Number of basis functions","Space (in Float64s)","Space to factorize")
        plotComparison(num_basis_fns[1:nmols],full_fact_error,"Full factorization",struct_fact_error,"Structured factorization","Number of basis functions","Error (1-norm)","Error to factorize")
        plotComparison(num_basis_fns[1:nmols],full_fact_fevals,"Full factorization",struct_fact_fevals,"Structured factorization","Number of basis functions","Fevals","Fevals to factorize")

        plotComparison(num_basis_fns[1:nmols],full_trans_time,"Full transformation",struct_trans_time,"Structured transformation","Number of basis functions","Time (seconds)","Time to transform")
        plotComparison(num_basis_fns[1:nmols],full_trans_space,"Full transformation",struct_trans_space,"Structured transformation","Number of basis functions","Space (in Float64s)","Space to transform")
        plotComparison(num_basis_fns[1:nmols],full_trans_error,"Full transformation",struct_trans_error,"Structured transformation","Number of basis functions","Error (1-norm)","Error to transform")

        plotComparison(num_basis_fns[1:nmols],full_fact_time+full_trans_time,"Full factorization plus transformation",struct_fact_time+struct_trans_time,"Structured factorization plus transformation","Number of basis functions","Time (seconds)","Time to factorize and transform")
    end

    if false
        figure()
        plot(1:maxn,ourtimes,color="red",label=l1)
        plot(1:maxn,comparedtimes,color="blue",label=l2)
        legend()
        xlabel("n=2 to 25")
        ylabel("Time in seconds")
        title("Time in seconds to tranform all integrals")
        savefig("timing.png")
        figure()
        plot(1:maxn,-log10(errors));
        xlabel("n=2 to 25")
        ylabel("-Log(error)")
        title("Frobenius norm error in transformed matrix")
        savefig("error.png")
        figure()
        plot(1:maxn,ranks);
        xlabel("n=2 to 25")
        ylabel("Numerical rank")
        title("Rank of the ERI matrix")
        savefig("rank.png")
        figure()
        plot(1:maxn,fevals);
        xlabel("n=2 to 25")
        ylabel("F evals")
        title("Number of two-electron integral evaluations")
        savefig("fevals.png")
    end
end

function plotComparison(x,y1,l1,y2,l2,xlab,ylab,titl)
    figure()
    plot(x,y1,"ro",label=l1)
    plot(x,y2,"bo",label=l2)
    legend()
    xlabel(xlab)
    ylabel(ylab)
    title(titl)
    savefig(string(titl,".png"))
end

function QsymTransform(L,piv,n)
    # Applies Qsym to G
    nsym,r = size(L)
    J = zeros(n^2,r)
    s = zeros(nsym)
    for i=1:r
        s[piv] = L[:,i]
        J[:,i] = structUnpack(s,n)
    end
    return J
end

function structUnpack(s,n)
    # Returns v=Qsym*s
    v = zeros(n^2)
    nsym = convert(Int64,n*(n+1)/2)
    for c=1:nsym
        i, j = unpackIndex(c,n)
        if i == j
            alpha = 1
            r = i+(j-1)*n
            v[r] = alpha*s[c]
        else
            alpha = 1/sqrt(2)
            r = i+(j-1)*n
            v[r] = alpha*s[c]
            r2 = j+(i-1)*n
            v[r2] = alpha*s[c]
        end
    end
    return v
end

function full_fact(A,tol)
    # Computes the full factorization of A
    full_fact_space = 0
    full_fact_fevals = 0
    n = size(A,1)
    col = zeros(n) # n vector
    full_fact_space += n
    L = Array{Float64,1}[]
    m = 1
    # evaluate diagonal
    d = diag(A) # n vector
    full_fact_space += n
    full_fact_fevals += n
    error = norm(d,1)
    piv = [1:n] # n vector
    full_fact_space += n
    while error > tol
        idx = indmax(d[piv[m:n]])
        i = idx+m-1
        piv[m], piv[i] = piv[i], piv[m]
        push!(L,zeros(n)) # n vector
        full_fact_space += n
        L[m][piv[m]] = sqrt(d[piv[m]])
        # evaluate column
        col = A[piv[1:n],piv[m]]
        full_fact_fevals += n
        for i=m+1:n
            s = 0
            for j=1:m-1
                s += L[j][piv[m]]*L[j][piv[i]]
            end
            L[m][piv[i]] = (col[i] - s)/L[m][piv[m]]
            d[piv[i]] = d[piv[i]] - L[m][piv[i]]*L[m][piv[i]]
        end
        error = sum(d[piv[m+1:n]])
        m += 1
    end
    L2 = zeros(n,m-1)
    full_fact_space += n*(m-1)
    for p=1:m-1
        L2[:,p] = L[p][piv]
    end
    return L2, piv, full_fact_space, full_fact_fevals
end

function full_trans(L,piv,C)
    # Computes the factorized full transformation of A
    n,n = size(C)
    n2,r = size(L)
    full_trans_space = 0
    B0 = zeros(n2,r)
    B0[piv,:] = L
    full_trans_space += n2*r
    B1 = C*reshape(B0,n,n*r)
    full_trans_space += n2*r
    B2 = C*reshape(permutedims(reshape(B1,n,n,r),[2 1 3]),n,n*r)
    full_trans_space += n2*r
    B3 = reshape(permutedims(reshape(B2,n,n,r),[2 1 3]),n^2,r)
    full_trans_space += n2*r
    return B3, full_trans_space
end

function struct_fact(A,tol)
    # Computes the structured factorization of A

    struct_fact_space = 0
    struct_fact_fevals = 0
    n2,n2 = size(A)
    n = convert(Int64,sqrt(n2))
    nsym = convert(Int64,n*(n+1)/2) # size(A,1)
    col = zeros(nsym) # nsym vector
    struct_fact_space += nsym
    L = Array{Float64,1}[]
    m = 1
    # evaluate diagonal
    # d = diag(A)
    d = struct_diag(A) # nsym vector
    struct_fact_space += nsym
    struct_fact_fevals += nsym
    error = norm(d,1)
    piv = [1:nsym] # nsym vector
    struct_fact_space += nsym
    while error > tol
        idx = indmax(d[piv[m:nsym]])
        i = idx+m-1
        piv[m], piv[i] = piv[i], piv[m]
        push!(L,zeros(nsym)) # n vector
        struct_fact_space += nsym
        L[m][piv[m]] = sqrt(d[piv[m]])
        # evaluate column
        # A[piv[1:nsym],piv[m]]
        col = struct_col(A,piv,m)
        struct_fact_fevals += nsym
        for i=m+1:nsym
            s = 0
            for j=1:m-1
                s += L[j][piv[m]]*L[j][piv[i]]
            end
            L[m][piv[i]] = (col[i] - s)/L[m][piv[m]]
            d[piv[i]] = d[piv[i]] - L[m][piv[i]]*L[m][piv[i]]
        end
        error = sum(d[piv[m+1:nsym]])
        m += 1
    end
    L2 = zeros(nsym,m-1)
    struct_fact_space += nsym*(m-1)
    for p=1:m-1
        L2[:,p] = L[p][piv]
    end
    return L2, piv, struct_fact_space, struct_fact_fevals
end

function unpackIndex(c,n)
    # Returns the indices (i,j) corresponding to packed index c
    # Assume that c = i-(1/2)(j-1)j+(j-1)n is between 1 and n(n+1)/2
    # where j is between 1 and n, and i is between j and n
    # Indexing from zero: c0 = i0-(1/2)j0(j0+1)+j0n, i0 >= j0
    # If i0 == j0, then c0 is quadratic in one variable k0 == i0 == j0
    # If i0 > j0 and we can find j0, all we need is the offset i0-j0 to find i0
    # Note: there is an interval of c0's with contant j0 (called interval-j0)
    # and i0 is determined by the offset within interval-j0
    c0 = c-1
    # Solve the quadratic equation k0^2-(2n+1)k0+2c0 = 0 for k0
    k0 = (1/2)*(2n+1-sqrt((2n+1)^2-8c0))
    # Interval-j0 is the integer component of k0
    j0 = floor(k0)
    # Let the first value in interval-j0 be called d0
    d0 = j0-(1/2)*j0*(j0+1)+j0*n
    # The offset i0-j0 is simply c0-d0
    i0 = j0+(c0-d0)
    return (i0+1,j0+1)
end

function struct_diag(A)
    # Evaluates the diagonal of Asym = Qsym'*A*Qsym
    n2,n2 = size(A)
    n = convert(Int64,sqrt(n2))
    nsym = convert(Int64,n*(n+1)/2)
    d = zeros(nsym)
    for c=1:nsym
        i, j = unpackIndex(c,n)
        r = i+(j-1)*n
        if i == j
            alpha = 1
        else
            alpha = 2
        end
        d[c] = alpha*A[r,r]
    end
    return d
end

function struct_col(A,piv,m)
    # Evaluates a column of Asym = Qsym'*A*Qsym
    n2,n2 = size(A)
    n = convert(Int64,sqrt(n2))
    nsym = convert(Int64,n*(n+1)/2)
    c = piv[m]
    i, j = unpackIndex(c,n)
    r = i+(j-1)*n
    # A[piv[1:n2],piv[r]]
    col = zeros(nsym)
    for p=1:nsym
        k, l = unpackIndex(piv[p],n)
        q = k+(l-1)*n
        if i == j && k == l
            alpha = 1
        elseif i == j && k > l || k == l && i > j
            alpha = sqrt(2)
        else
            alpha = 2
        end
        col[p] = alpha*A[q,r]
    end
    return col
end

function struct_trans(L,piv,C)
    # Computes the factorized structured transformation of A
    n,n = size(C)
    struct_trans_space = 0
    J = QsymTransform(L,piv,n)
    n2,r = size(J)
    struct_trans_space += n2*r
    B1 = C*reshape(J,n,n*r)
    struct_trans_space += n2*r
    B2 = C*reshape(permutedims(reshape(B1,n,n,r),[2 1 3]),n,n*r)
    struct_trans_space += n2*r
    B3 = reshape(permutedims(reshape(B2,n,n,r),[2 1 3]),n^2,r)
    struct_trans_space += n2*r
    return B3, struct_trans_space
end

function CompareIntegralTransformation(n)
    # Compares SVDIntegralTransformation with FastIntegralTransformation
    println("n= ",n)

    # Random n^2 by n^2 ERI matrix of rank 2n with correct symmetries
    # Unfolded as A_ij,kl=(ij|kl)
    r = 2n
    A = RandERIMatrix(n,r);
    r = rank(A)
    # Coefficient matrix from Hartree-Fock method
    C = randn(n,n);
    # tol determines when to stop LazyLDLT
    tol = 1e-13;
    # Reshapes ERI tensor as T[i,j,k,l]=(ij|kl)
    T = reshape(A,n,n,n,n)

    tic()
    M1 = StateOfTheArtIntegralTransformation(T,C)
    println("SOTA total time")
    naivetime = toc();

    tic()
    M2,sig,Z = FastIntegralTransformation(A,C,tol);
    println("FIT total time")
    fasttime = toc();

    error = normfro(M1-M2);
    println("Error: ",error)
    return naivetime, fasttime, error, r
end

function RandERIMatrix(n,r)
    # Returns a random positive semidefinite matrix unfolded from the ERI
    # T(i1,i2,i3,i4) as A(i1+(i2-1)n,i3+(i4-1)n)
    A = zeros(n^2,n^2)
    for i=1:r
        B = randn(n,n)
        B = (B + B')/2
        A += vec(B)*vec(B)'
    end
    return A
end

function StateOfTheArtIntegralTransformation(T,C)
    # Returns the transformed integrals by a series of modal products
    n,n,n,n = size(T)
    
    println("SOTA step 1 time")
    tic()
    temp1 = zeros(n,n^3)
    temp1 = C'*reshape(T,n,n^3)
    toc()
    println("SOTA step 2 time")
    tic()
    temp2 = zeros(n,n^3)
    temp2 = C'*reshape(permutedims(reshape(temp1,n,n,n,n),[2 1 3 4]),n,n^3)
    toc()
    println("SOTA step 3 time")
    tic()
    temp3 = zeros(n,n^3)
    temp3 = C'*reshape(permutedims(reshape(temp2,n,n,n,n),[3 2 1 4]),n,n^3)
    toc()
    println("SOTA step 4 time")
    tic()
    temp4 = zeros(n,n^3)
    temp4 = C'*reshape(permutedims(reshape(temp3,n,n,n,n),[4 2 3 1]),n,n^3)
    MO = zeros(n^2,n^2)
    MO = reshape(permutedims(reshape(temp4,n,n,n,n),[2 3 4 1]),n^2,n^2)
    toc()
    return MO
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
                                    M[p,q,r,s] += C[mu,p]*C[nu,q]*C[lam,r]*C[sig,s]*T[mu,nu,lam,sig]
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
    n = convert(Int64,sqrt(n2))
    # Computes the Kronecker Product SVD using lazy-evaluation LDL^T
    tic()
    (B, sig) = StructLDLT(A,tol);
    println("Factorization time")
    toc()
    tic()
    rsym = length(sig)
    #Z = cell(rsym)
    Z = zeros(n,n)
    M = zeros(n2,n2)
    #Mt = zeros(n2,n2)
    for k=1:rsym
        #Z[k] = (C'*(B[k]*C))
        Z = (C'*(B[k]*C))
        #v = vec(Z[k])
        #M = M + sig[k]*v*v'
    end
    println("Transformation time")
    toc()
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
    println("LazyLDLT time")
    Asymtemp = deepcopy(Asym)
    tic()
    Lsym, psym, r, _ = LAPACK.pstrf!('L',Asymtemp,1e-6)
    #    (Lsym,sig,psym) = LazyLDLT(F, tol, nsym, Asym);
    toc()
    sig = diag(Lsym)
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
    #global numfevals = 0;
    v = zeros(n);
    w = zeros(n);
    temp = zeros(n);
    d = zeros(n);
    L = eye(n,n);
    piv = [i for i in 1:n];
    thresh = n*d[1]*tol
    # Get the diagonal elements
    d = diag(A)
    #for i=1:n
    #    d[i] = F(i,i,A);
    #end
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
        w[j+1:n] = A[piv[j+1:n],piv[j]]
        #    temp[j+1:n] = A[piv[j+1:n],piv[j]]
        #    println("Length of w: ", length(w))
        #    println("Length of piv[j+1:n]: ", length(piv[j+1:n]))
        #    println("Length of j+1:n: ", length(j+1:n))
        #for i=j+1:n
        #    w[i] = F(piv[i],piv[j],A);
        #end
        # Compute d(j) and L(:,j)
        v[1:j-1] = L[j,1:j-1].*d[1:j-1]'
        #    temp[1:j-1] = L[j,1:j-1].*d[1:j-1]'
        #    println("Length of v: ", length(v))
        #    println("Length of 1:j-1: ", length(1:j-1))
        #    for i=1:j-1
        #        v[i]= L[j,i]*d[i];
        #    end
        #    println("Error: ",norm(temp[1:j-1]-v[1:j-1]))
        c = L[j,1:j-1]*v[1:j-1]
        d[j] = alpha - c[1]
        L[j+1:n,j] = (w[j+1:n] - L[j+1:n,1:j-1]*v[1:j-1])/d[j]
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

function TestCholeskyTransformation(A,C)
    # Returns the transformed integrals by a series of modal products
    n2,r = size(A)
    n = convert(Int64,sqrt(n2))
    #    temp1 = zeros(n,n*r)
    temp1 = C'*reshape(A,n,n*r)
    #    temp2 = zeros(n,n*r)
    temp2 = C'*reshape(permutedims(reshape(temp1,n,n,r),[2 1 3]),n,n*r)
    MO = zeros(n^2,r)
    MO = reshape(permutedims(reshape(temp2,n,n,r),[2 1 3]),n^2,r)
    return MO
end

function compute_diag(A)
    return deepcopy(diag(A))
end

function compute_col(A,j)
    #    println("j: ",j)
    #    println("nsym: ",nsym)
    #    println("size of A: ",size(A))
    return deepcopy(A[:,j])
end

function cholesky2(A,n,delta)
    # initial dimensions
    #nsym = convert(Int64,nb*(nb+1)/2)
    q = 1
    # compute_diagonal
    d = compute_diag(A)
    #    display(d)
    # temporary cholesky factor
    #    L = zeros(nsym,nsym)
    L = Array{Float64,1}[]
    #    println("Size of L: ",size(L))
    # list of selected pivots
    #    pivots = 1:nsym
    pivots = Int64[]
    # cholesky procedure
    piv = [i for i in 1:n];
    thresh = n*d[1]*tol
    # Get the diagonal elements
    #d = diag(A)
    #for i=1:n
    #    d[i] = F(i,i,A);
    #end
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
        #    L[k,:], L[j,:] = L[j,:], L[k,:];
        # L(:,[k j]) = L(:,[j k]);
        #    L[:,k], L[:,j] = L[:,j], L[:,k];
        #    while q <= nb
        #        println("q: ",q)
        #  select the pivot
        #        dvec = d[q:nb]
        #        println("d: ",d)
        #        (dmax,idx) = maximum(dvec), indmax(dvec)
        #        pivot = idx+q-1
        #        pivot = 1
        #        dmax = d[1]
        #        for p=1:nb
        #            if dmax < d[p]
        #                dmax = d[p]
        #                pivot = p
        #            end
        #        end

        #pivot = 1
        #dmax = d[1]
        #for p=1:nsym
        #    if dmax < d[p]
        #        dmax = d[p]
        #        pivot = p
        #    end
        #end
        # Swap dk and dj
        # d([k j]) = d([j k]);
        #        d[q], d[pivot] = d[pivot], d[q]
        # Update piv
        # piv([k j]) = piv([j k]);
        #        pivots[pivot], pivots[p] = pivots[q], pivots[pivot]

        #  check to see if convergence reached
        #        if dmax < delta || dmax < 0
        #            println("dmax: ",dmax)
        #            println("delta: ",delta)
        #            println("Convergence reached")
        #            break
        #        end
        #  if here, we're trying to add this row
        #        pivots[q] = pivot
        #        push!(pivots,pivot)
        #        println("Pivot added: ",pivot)
        #        println("Pivots: ",pivots)
        #        Lqq = sqrt(dmax)
        #  if here, we're really going to add this row
        #  (m|Q)
        #        L[:,q] = compute_col(pivot,nsym)
        push!(L,compute_col(A,pivot))

        #        println("Size of Lq: ",size(L[q]))
        #        println("Lq: ",L[q])
        #  [(m|Q) - L_m^P L_Q^P]
        for p=1:q
            L[q] = L[q] - L[p][pivots[q]]*L[p]
        end
        #  1/L_QQ [(m|Q) - L_m^P L_Q^P]
        L[q] = L[q]/Lqq
        #  zero the upper triangle
        for p=1:q
            L[q][pivots[p]] = 0
        end
        # set the pivot factor
        L[q][pivot] = Lqq
        #  update the schur complement diagonal
        for p=1:nb
            d[p] -= L[q][p]*L[q][p]
        end
        #  force truly zero elements to be zero
        for p=1:q
            d[pivots[p]] = 0
        end
        q += 1
        #    c = L[j,1:j-1]*v[1:j-1]
        #    d[j] = alpha - c[1]
        #    L[j+1:n,j] = (w[j+1:n] - L[j+1:n,1:j-1]*v[1:j-1])/d[j]
        #    crit = d[j]
    end
    # copy into a more permament matrix object
    println("q: ",q)
    G = zeros(nb,q-1)
    for p=1:q-1
        G[:,p] = L[p]
    end
    return G, pivots
end

function inv_col(k,n)
    # function i = inv_col(k,n)
    # n is a d-vector with positive integer entries
    # k is an integer that satisfies 1<=k<=prod(n)
    # i is a d-vector with positive integer entries so
    # k = i(1) + (i(2)-1)n(1) + (i(3)-1)n(1)n(2) + ... + (i(d)-1)n(1)...n(d-1)
    # GVL4: Section 12.4.5
    d = length(n);
    i = zeros(d,1);
    N = prod(n);
    k = k-1;
    for j=d:-1:1
        N = N/n[j];
        f = floor(k/N);
        k = k-N*f;
        i[j] = f+1;   
    end
    return i
end

function unpack(AP,n)
    # Converts AP in packed storage format to A in full format
    A = zeros(n,n)
    for j=1:n
        k = j*(j-1)/2
        A[1:j,j] = AP[1+k:j+k]
        A[j,1:j-1] = AP[1+k:j-1+k]'
    end
    return A
end

function pack(A,n)
    # Converts full A to packed storage
    v = zeros(convert(Int64,n*(n+1)/2))
    for j=1:n
        k = j*(j-1)/2
        v[1+k:j+k] = A[1:j,j]
    end
    return v
end

function LazyCholesky(A,tol)
    # Source: "On the low rank approximation by the pivoted cholesky" by Harbrecht, Peters, and Schneider
    # Devectorized
    n = size(A,1)
    col = zeros(n) # n vector
    L = Array{Float64,1}[]
    m = 1;
    # evaluate diagonal
    d = diag(A) # n vector
    error = norm(d,1)
    piv = [1:n] # n vector
    while error > tol
        idx = indmax(d[piv[m:n]])
        i = idx+m-1
        piv[m], piv[i] = piv[i], piv[m]
        push!(L,zeros(n)) # n vector
        L[m][piv[m]] = sqrt(d[piv[m]])
        # evaluate column
        col = A[piv[1:n],piv[m]]
        for i=m+1:n
            s = 0
            for j=1:m-1
                s += L[j][piv[m]]*L[j][piv[i]]
            end
            L[m][piv[i]] = (col[i] - s)/L[m][piv[m]]
            d[piv[i]] = d[piv[i]] - L[m][piv[i]]*L[m][piv[i]]
        end
        error = sum(d[piv[m+1:n]])
        m += 1
    end
    G = zeros(n,m-1)
    for p=1:m-1
        G[:,p] = L[p][piv]
    end
    return (G,piv)
end

function GaxpyLazyCholesky(A,tol)
    # Source: https://github.com/psi4/psi4public/blob/c98c9d3f28c6899c6cd8ed3092ed6887c2426bd9/src/lib/lib3index/cholesky.cc
    # Vectorized
    n = size(A,1)
    L = Array{Float64,1}[]
    m = 1;
    # evaluate diagonal
    d = diag(A)
    error = norm(d,1)
    piv = [1:n]
    col = zeros(n)
    while error > tol
        idx = indmax(d[piv[m:n]])
        i = idx+m-1
        piv[m], piv[i] = piv[i], piv[m]
        push!(L,zeros(n))
        L[m][piv[m]] = sqrt(d[piv[m]])
        # evaluate column
        col[m+1:n] = A[piv[m+1:n],piv[m]]
        L[m][piv[m+1:n]] = col[m+1:n]
        for i=1:m-1
            L[m][piv[m+1:n]] = L[m][piv[m+1:n]] - L[i][piv[m]]*L[i][piv[m+1:n]]
        end
        L[m][piv[m+1:n]] /= L[m][piv[m]]
        d[piv[m+1:n]] = d[piv[m+1:n]] - L[m][piv[m+1:n]].*L[m][piv[m+1:n]]
        error = sum(d[piv[m+1:n]])
        m += 1
    end
    G = zeros(n,m-1)
    for p=1:m-1
        G[:,p] = L[p][piv]
    end
    return (G,piv)
end

function Qnn(n)
    # Q is the n^2 by n^2 orthgonal matrix which 
    # block diagonalizes symmetric perfect shuffle invariant matrices

    Q = zeros(n^2,n^2)
    e = eye(n,n)
    k = 0
    # Define Qsym
    for i=1:n
        for j=i:n
            k = k+1
            if i == j
                Q[:,k] = kron(e[:,i],e[:,i])
            else
                Q[:,k] = (kron(e[:,i],e[:,j]) + kron(e[:,j],e[:,i]))/sqrt(2)
            end
        end
    end
    # Define Qskew
    for i=1:n
        for j=i+1:n
            k = k+1
            Q[:,k] = (kron(e[:,i],e[:,j]) - kron(e[:,j],e[:,i]))/sqrt(2)
        end
    end
    return Q
end

function start()

    TestIntegralTransformation()
    #n = 60
    #C = randn(n,n)
    #r = 6*n
    #QPL = randn(n^2,r)
    #tic()
    #test1 = TestCholeskyTransformation(QPL,C)
    #toc()
    #tic()
    #test2 = kron(C',C')*QPL
    #toc()
    #error = norm(test1-test2)
    #n = 3
    #A = RandERIMatrix(n,n)
    #A = randn(3,3)
    #A = A'*A
    #Ain = deepcopy(A)
    #G,piv = LazyCholesky(Ain,1e-6)
    #L = zeros(n,n)
    #L[piv,:] = G
    #G[piv,:] = G
    #error = norm(A - G*G',1)
    #println("Error: ",error)
    #n = 10
    if false
        #v = randn(convert(Int64,n*(n+1)/2))
        #display(unpack(v,n))
        #println()
        #error = norm(v-pack(unpack(v,n),n))
        #println("Pack error: ",error)
        #nsym = convert(Int64,n*(n+1)/2)
        n = 100
        #global A
        A = randn(n-3,n)
        A = A'*A
        println("Rank of A: ",rank(A))
        #A = randn(n,n-3)
        #println(size(A))
        #A = A*A'
        #println(size(A))
        println("Size of A: ",size(A))
        Ain = deepcopy(A)
        tic()
        G,piv = GaxpyLazyCholesky(Ain,1e-6)
        toc()
        #println("G: ")
        #display(G)
        #println()
        #println("piv: ",piv)
        Gaxpyerror = norm(A[piv,piv] - G*G')
        Ain = deepcopy(A)
        tic()
        G,piv = LazyCholesky(Ain,1e-6)
        toc()
        #println("G: ")
        #display(G)
        #println()
        #println("piv: ",piv)
        Devecerror = norm(A[piv,piv] - G*G')
        println("Gaxpy error: ",Gaxpyerror,"\nDevec error: ",Devecerror)
        #println(piv)
        #println(size(G))
        #println(size(piv))
        #Ain = deepcopy(A)
        #L, d, piv = LazyLDLT(F,1e-6,n,Ain)
        #error = norm(A[piv,piv] - L*diagm(d)*L')
        #error = norm(A[piv,piv] - G*G')
        #println("Error: ",error)
        n = 3
        A = randn(n^2,n^2)
        C = randn(n,n)
        A = A'*A
        Ain = deepcopy(A)
        G,piv = LazyCholesky(A,1e-6)
        L = zeros(n^2,n^2)
        L[piv,:] = G
        B,_ = full_trans(L,C)
        error = norm(kron(C,C)*A*kron(C,C)' - B*B')
        println("Error: ",error)
    end

    n = 3
    r = 3
    Q = Qnn(n)
    A = randn(9,9)
    A = A + A[PerfShuff(3,3),PerfShuff(3,3)]
    println("A: perfect shuffle invariant")
    display(A)
    println("\nQ'AQ: ")
    display(Q'*A*Q)
    B = randn(9,9)
    B = B + B'
    B = B + B[PerfShuff(3,3),PerfShuff(3,3)]
    println("\nB: symmetric perfect shuffle invariant")
    display(B)
    println("\nQ'BQ: ")
    display(Q'*B*Q)
    C = RandERIMatrix(n,r)
    println("\nC: [12]x[34] unfolding of ((ij)(kl))-symmetric tensor")
    display(C)
    println("\nQ'CQ: ")
    display(Q'*C*Q)
    println("\nRank of Q'CQ: ",rank(Q'*C*Q))
    P = eye(n^2,n^2)
    P = P[:,PerfShuff(n,n)]
    D = kron(eye(n,n),P,eye(n,n))*vec(C)
    D = reshape(D,n^2,n^2)
    println("\nD: [13]x[24] unfolding of ((ij)(kl))-symmetric tensor")
    display(D)
    println("\nQ'DQ: ")
    display(Q'*D*Q)
    println("\nRank of Q'DQ: ",rank(Q'*D*Q))
    println()
    #L, piv, _, _ = struct_fact(A,1e-6)
    #J = QsymTransform(L,piv,n)
    #error = norm(A-J*J')
    #println("Error: ",error)
    #println("A: ")
    #display(A)
    #println("\nJJ': ")
    #display(J*J')
    #println("Q: ")
    #Q = Qnn(n)
    #nsym = convert(Int64,n*(n+1)/2)
    #Qplus = Q[:,1:nsym]
    #Aplus = Qplus'*A*Qplus
    #println("\nAplus: ")
    #display(Aplus)
    #error = norm(Aplus[piv,piv]-L*L')
    #Aplus[piv,piv] = L*L'
    #println("\n(P'L)(P'L)': ")
    #display(Aplus)
    #println("\nStructured factorization error: ",error)

    #s = randn(nsym)
    #v = structUnpack(s,n)
    #println("v: ",v)
    #w = Qplus*s
    #println("w: ",w)
    #error = norm(w-v)
    #println("Error: ",error)

end

TestIntegralTransformation()