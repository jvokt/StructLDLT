function SizeRankVersusSpeedup()
    # table of n = 500, 1000, 1500, 2000; r = n, n/2, n/10; f = 0, 1
    # ratio of L = UnStructPerfShuff(A) time divided by L1,L2 = StructPerfShuff(A) time
    
    tol = 1e-6
    trials = 50
    n_range = [20, 40, 60, 80]#, 100]
    unstruct_times = zeros(length(n_range),5,trials)
    struct_times = zeros(length(n_range),5,trials)
    for n_i = 1:length(n_range)
        n = n_range[n_i]
        nsym = convert(Int64,n*(n+1)/2)
        nskew = convert(Int64,n*(n-1)/2)
        Q = sparse(Qnn(n))
        Qsym = Q[:,1:nsym]
        Qskew = Q[:,nsym+1:end]
        r_range = [(nsym,nskew), (nsym,n/2), (nsym,0), (n/2,n/2), (n,0)]
        for r_i = 1:length(r_range)
            r1,r2 = convert((Int64,Int64),r_range[r_i])
            r = r1+r2
            println("n: ",n,", r1: ",r1,", r2: ",r2)
            B,C = RandPerfShuff(n,r1,r2)
            # A = FullPerfShuff(B,C)
            X = Qsym*B*Qsym'
            Y = Qskew*C*Qskew'
            A = X + Y
            for t_i = 1:trials
                Atemp = deepcopy(A)
                Btemp = deepcopy(B)
                Ctemp = deepcopy(C)

                tic() # <- start
                L,piv,rank,_ = UnStructPerfShuff(Atemp)
                unstruct_times[n_i,r_i,t_i] = toc() # <- stop
                # L[piv,:] = tril(L)
                # L = L[:,1:rank]
                # println("error: ",norm(A-L*L'))
                
                tic() # <- start
                Lsym,psym,rsym,Lskew,pskew,rskew = StructPerfShuff(Btemp,Ctemp)
                struct_times[n_i,r_i,t_i] = toc() # <- stop
                # Lsym[psym,:] = tril(Lsym)
                # Gsym = Lsym[:,1:rsym]
                # Gsym = Qsym*Gsym
                # Lskew[pskew,:] = tril(Lskew)
                # Gskew = Lskew[:,1:rskew]
                # Gskew = Qskew*Gskew
                # println("error: ",norm(A-Gsym*Gsym'-Gskew*Gskew'))
            end
        end
    end
    display(mean(unstruct_times,3)./mean(struct_times,3))
    println()
    display(mean(unstruct_times./struct_times,3))
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

function Onn(n)
    # O is the n^2 by n^2 zero-one matrix which 
    # block diagonalizes symmetric perfect shuffle invariant matrices

    O = zeros(n^2,n^2)
    e = eye(n,n)
    k = 0
    # Define Osym
    for i=1:n
        for j=i:n
            k = k+1
            if i == j
                O[:,k] = kron(e[:,i],e[:,i])
            else
                O[:,k] = (kron(e[:,i],e[:,j]) + kron(e[:,j],e[:,i]))
            end
        end
    end
    # Define Oskew
    for i=1:n
        for j=i+1:n
            k = k+1
            O[:,k] = (kron(e[:,i],e[:,j]) - kron(e[:,j],e[:,i]))
        end
    end
    return O
end

function tnn(n)
    unity = [j+(j-1)*n-(j-1)j/2 for j=1:n]
    t = ones(n^2)/sqrt(2)
    t[unity] = 1
    return t
end

function SymUnpackRows(B,n)
    # [unpackIndex(c,n) for c=1:n*(n+1)/2]
    G = zeros(n^2,size(B,2))
    for c=1:n*(n+1)/2
        i, j = unpackIndex(c,n)
        if i > j
            r1 = i+(j-1)*n
            G[r1,:] = B[c,:]/sqrt(2)
            r2 = j+(i-1)*n
            G[r2,:] = B[c,:]/sqrt(2)
        else
            r1 = i+(j-1)*n
            G[r1,:] = B[c,:]
        end            
    end
    return G
end

function CompareSymUnpackRowsQsym()
    for n=10:10:100
        nsym = convert(Int64,n*(n+1)/2)
        B = randn(nsym,n)
        Q = sparse(Qnn(n))
        Qsym = Q[:,1:nsym]
        
        tic()
        G1 = Qsym*B
        toc()

        tic()
        G2 = SymUnpackRows(B,n)
        toc()
        
        println("Error: ",norm(G1-G2))
    end
end

function RandPerfShuff(n,r1,r2)
    # Returns a rank-r matrix with Perfect Shuffle Symmetry
    nsym = convert(Int64,n*(n+1)/2)
    nskew = convert(Int64,n*(n-1)/2)
    B = randn(nsym,r1)
    B = B*B'
    C = randn(nskew,r2)
    C = C*C'
    return B,C
end

function lapack_chol(A,tol=1e-6)
    A, piv, rank, info = LAPACK.pstrf!('L', A, tol)
    A[piv,:] = tril(A)
    return A[:,1:rank]
end

function FullPerfShuff(B,C)
    Q = sparse(Qnn(n))
    Qsym = Q[:,1:nsym]
    Qskew = Q[:,nsym+1:end]
    X = Qsym*B*Qsym'
    Y = Qskew*C*Qskew'
    return X + Y
end

function UnStructPerfShuff(A,tol=1e-6)
    # Computes the rank-revealing Cholesky factorization of A, without utilizing Perfect Shuffle Symmetry Structure
    return LAPACK.pstrf!('L', A, tol)
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

function struct_block(A,lower_i,lower_j,block_size_i,block_size_j)
    # Evaluates a column of Asym = Qsym'*A*Qsym
    n2,n2 = size(A)
    n = convert(Int64,sqrt(n2))
    nsym = convert(Int64,n*(n+1)/2)
    block = zeros(block_size_i,block_size_j)
    # A[lower_i:upper_i,lower_j:upper_j]
    for c=1:block_size_j
        row = c + lower_j - 1
        i, j = unpackIndex(row,n)
        r = i+(j-1)*n
        for p=1:block_size_i
            col = p + lower_i - 1
            k, l = unpackIndex(col,n)
            q = k+(l-1)*n
            if i == j && k == l
                alpha = 1
            elseif i == j && k > l || k == l && i > j
                alpha = sqrt(2)
            else
                alpha = 2
            end
            block[p,c] = alpha*A[q,r]
        end
    end
    return block
end

function StructPerfShuff(B,C,tol=1e-6)
    # Computes the rank-revealing Cholesky factorization of A, utilizing PerfShuffSymmetric structure
    Asym,psym,rsym,_ = LAPACK.pstrf!('L', B, tol)
    Askew,pskew,rskew,_ = LAPACK.pstrf!('L', C, tol)
    return Asym,psym,rsym,Askew,pskew,rskew
end

function PerfShuff(p,r)
# function v = PerfShuff(p,r)
# Vector representation of the perfect shuffle permutation.
# GVL4: Setion 1.2.11
# p and r are positive integers. If n = pr then v is a length-n row
#   vector with the property that if x is a column n-vector and
#   y = x(v), then y = I(v,:)*x = P_{p,r}* x is the mod-p perfect shuffle
#   of x.
n = p*r;
v = zeros(n);
for k=1:r
    v[(k-1)*p+1:k*p] = k:r:n;
end
return v
end

function TestRank()
    f = 0
    n_range = [10]#, 1000, 1500, 2000]
    for n_i = 1:length(n_range)
        n = n_range[n_i]
        r_range = [n^2, n^2/2, n^2/10]
        for r_i = 1:length(r_range)
            r = r_range[r_i]
            println("n: ",n,", r: ",r)
            A = RandPerfShuffSym(n,r,f)
            nsym = convert(Int64,n*(n+1)/2)
            B11 = struct_block(A,1,1,nsym,nsym)
            println("rank(B11): ", rank(B11))
            B22 = struct_block(A,1,1,nsym,nsym)
            println("rank(B22): ",rank(B22))
        end
    end
end

function TestPerfShuffSymmetries()
    P = eye(9,9)[:,PerfShuff(3,3)]
    A = randn(n,n)
    A = A + A'
    A = A + P*A*P
    T_A = reshape(A,3,3,3,3)
    B = RandPerfShuffSym(3,6,0)
    T_B = reshape(B,3,3,3,3)
    for p in permutations([1 2 3 4])
        println(norm(vec(T_A-permutedims(T_A,p))))
    end
    
end

SizeRankVersusSpeedup()
#TestRank()
#CompareSymUnpackRowsQsym()