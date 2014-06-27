function SizeRankVersusSpeedup()
    # table of n = 500, 1000, 1500, 2000; r = n, n/2, n/10; f = 0, 1
    # ratio of L = UnStructCentro(A) time divided by L1,L2 = StructCentro(A) time

    tol = 1e-6
    trials = 50
    n_range = [500, 1000, 1500, 2000]
    unstruct_times = zeros(length(n_range),5,trials)
    struct_times = zeros(length(n_range),5,trials)
    for n_i = 1:length(n_range)
        n = n_range[n_i]
        r_range = [(n/2,n/2), (n/2,n/10), (n/2,0), (n/100,n/100), (n/100,0)]
        for r_i = 1:length(r_range)
            r1,r2 = r_range[r_i]
            println("n: ",n,", r1: ",r1,", r2: ",r2)
            A = RandCentro(n,r1,r2)
            L = deepcopy(A)
            m = convert(Int64,n/2)
            A11 = A[1:m,1:m]
            A12 = A[1:m,m+1:n]
            B1 = A11 + A12[:,m:-1:1]
            B2 = A11 - A12[:,m:-1:1]
            for t_i = 1:trials
                tic() # <- start
                L,piv,rank,_ = LAPACK.pstrf!('L', L, tol)
                unstruct_times[n_i,r_i,t_i] = toc() # <- stop
                #L[piv,:] = tril(L)
                #G = L[:,1:rank]
                # println("error: ",norm(A-G*G'))
                # L = UnStructCentro(Acopy,tol)
                
                tic() # <- start
                B1,piv1,rank1,_ = LAPACK.pstrf!('L', B1, tol)
                B2,piv2,rank2,_ = LAPACK.pstrf!('L', B2, tol)
                struct_times[n_i,r_i,t_i] = toc() # <- stop
                #B1[piv1,:] = tril(B1)
                #B2[piv2,:] = tril(B2)
                #G1 = [B1[:,1:rank1]; B1[m:-1:1,1:rank1]]/sqrt(2)
                #G2 = [B2[:,1:rank2]; -B2[m:-1:1,1:rank2]]/sqrt(2)
                #println("error: ",norm(A-G1*G1'-G2*G2'))
                # L1 = UnStructPerfShuff(Asym,tol)
                # L2 = UnStructPerfShuff(Askew,tol)
            end
        end
    end
    display(mean(unstruct_times,3)./mean(struct_times,3))
    println()
    display(mean(unstruct_times./struct_times,3))
end

function RandCentro(n,r1,r2)
    # Returns a rank-r matrix with Centrosymmetry
    A = zeros(n,n)
    for i=1:r1
        v = randn(n,1)
        v += v[n:-1:1]
        A += v*v'
    end
    for i=1:r2
        v = randn(n,1)
        v -= v[n:-1:1]
        A += v*v'
    end
    return A
end

function lapack_chol(A,tol)
    Aout, piv, rank, info = LAPACK.pstrf!('L', A, tol)
    L = tril(Aout)
    L[piv,1:rank] = L[:,1:rank]
    return L[:,1:rank]
end

function UnStructCentro(A,tol)
    Aout, piv, rank, info = LAPACK.pstrf!('L', A, tol)
    L = tril(Aout)
    L[piv,1:rank] = L[:,1:rank]
    return L[:,1:rank]
end

function UnStructCentro2(A,tol,block_size)
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

function StructCentro(A,tol)
    # Computes the rank-revealing Cholesky factorization of A, utilizing Centrosymmetric structure
    n = size(A,1)
    m = convert(Int64,n/2)
    A11 = A[1:m,1:m]
    A12 = A[1:m,m+1:n]
    B11 = A11 + A12[:,m:-1:1]
    L1 = UnStructCentro(A11,tol)
    B22 = A11 - A12[:,m:-1:1]
    L2 = UnStructCentro(A12,tol)
    return L1,L2
end

function TestRank()
    ranks = zeros(4,3,2,2)
    n_range = [500]#, 1000, 1500, 2000]
    for n_i = 1:length(n_range)
        n = n_range[n_i]
        r_range = [n, n/2, n/10]
        for r_i = 1:length(r_range)
            r = r_range[r_i]
            f_range = [0,1]
            for f_i = 1:length(f_range)
                f = f_range[f_i]
                println("n: ",n,", r: ",r,", f: ",f)
                A = RandCentro(n,r,f)
                m = convert(Int64,n/2)
                A11 = A[1:m,1:m]
                A12 = A[1:m,m+1:n]
                B11 = A11 + A12[:,m:-1:1]
                println("rank(B11): ", rank(B11))
                B22 = A11 - A12[:,m:-1:1]
                println("rank(B22): ",rank(B22))
            end
        end
    end
end

SizeRankVersusSpeedup()
#TestRank()