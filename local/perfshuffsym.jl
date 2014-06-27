function SizeRankVersusSpeedup()
    # table of n = 500, 1000, 1500, 2000; r = n, n/2, n/10; f = 0, 1
    # ratio of [L,D,P] = UnStructPerfShuff(A) divided by [Lp,Dp,Pp,Lm,Dm,Pm] = StructPerfShuff(A)

    tol = 1e-6

    f = 0
    n_range = [10, 20, 30, 40]
    unstruct_times = zeros(length(n_range),3)
    struct_times = zeros(length(n_range),3)
    for n_i = 1:length(n_range)
        n = n_range[n_i]
        r_range = [n^2, n^2/2, n^2/10]
        for r_i = 1:length(r_range)
            r = r_range[r_i]
            println("n: ",n,", r: ",r)
            A = RandPerfShuffSym(n,r,f)
            tic()
            L = UnStructPerfShuff(A,tol,192)
            unstruct_times[n_i,r_i] = toc()
            tic()
            Lp,Lm = StructPerfShuff(A,tol)
            struct_times[n_i,r_i] = toc()
        end
    end
    display(unstruct_times./struct_times)
end

function RandPerfShuffSym(n,r,f)
    # Returns a rank-r matrix with Perfect Shuffle Symmetry
    A = zeros(n^2,n^2)
    for i=1:r
        B = randn(n,n)
        if f == 1 && i > n*(n+1)/2
            B -= B'
        else
            B += B'
        end
        A += vec(B)*vec(B)'
    end
    return A
end

function lapack_chol(A,tol)
    Aout, piv, rank, info = LAPACK.pstrf!('L', A, tol)
    L = tril(Aout)
    L[piv,1:rank] = L[:,1:rank]
    return L[:,1:rank]
end

function UnStructPerfShuff(A,tol,block_size)
    Aout, piv, rank, info = LAPACK.pstrf!('L', A, tol)
    L = tril(Aout)
    L[piv,1:rank] = L[:,1:rank]
    return L[:,1:rank]
end

function UnStructPerfShuff2(A,tol,block_size)
    n = size(A,1)
    num_blocks = convert(Int64,div(n,block_size)) + (convert(Bool,mod(n,block_size)) ? 1 : 0)
    D = Array{Float64,2}[]
    for i=1:num_blocks
        range_i = (i-1)*block_size+1:min(i*block_size,n)
        push!(D,A[range_i,range_i])
    end
    L = [Array{Float64,2}[] for i=1:num_blocks]
    piv = [1:num_blocks]
    error = sum(map(trace,D))
    j = 1
    while error > tol
        idx = indmax(map(trace,D[piv[j:num_blocks]]))+j-1
        piv[j],piv[idx] = piv[idx],piv[j]
        lower_j = (piv[j]-1)*block_size+1
        upper_j = min(piv[j]*block_size,n)
        range_j = lower_j:upper_j
        block_size_j = upper_j-lower_j+1
        G = lapack_chol(D[piv[j]],tol)
        r = size(G,2)
        # L[piv[j]][j] = G
        push!(L[piv[j]],G)
        for i=j+1:num_blocks
            lower_i = (piv[i]-1)*block_size+1
            upper_i = min(piv[i]*block_size,n)
            range_i = lower_i:upper_i
            block_size_i = upper_i-lower_i+1
            # GEMM update of subdiagonal blocks
            S = A[range_i,range_j]
            for k=1:j-1
                S -= L[piv[i]][k]*L[piv[j]][k]'
            end
            G = S/L[piv[j]][j]'
            # SYRK update of diagonal blocks
            D[piv[i]] -= G*G'
            # L[piv[i]][j] = G
            push!(L[piv[i]],G)
        end
        error = sum(map(trace,D[piv[j+1:num_blocks]]))
        j+=1
    end
    L2 = zeros(n,min(block_size*(j-1),n))
    for jj=1:j-1
        lower_j = (jj-1)*block_size+1
        for ii=jj:num_blocks
            lower_i = (piv[ii]-1)*block_size+1
            block_size_i,block_size_j = size(L[piv[ii]][jj])
            range_i = lower_i:min(lower_i+block_size_i-1,n)
            range_j = lower_j:min(lower_j+block_size_j-1,n)
            L2[range_i,range_j] = L[piv[ii]][jj]
        end
    end
    return L2
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

function StructPerfShuff(A,tol)
    # Computes the rank-revealing Cholesky factorization of A, utilizing PerfShuffsymmetric structure
    n2 = size(A,1)
    n = convert(Int64,sqrt(n2))
    nsym = convert(Int64,n*(n+1)/2)
    B11 = struct_block(A,1,1,nsym,nsym)
    L1 = UnStructPerfShuff(B11,tol,192)
    B22 = struct_block(A,1,1,nsym,nsym)
    L2 = UnStructPerfShuff(B22,tol,192)
    return L1,L2
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

SizeRankVersusSpeedup()
#TestRank()