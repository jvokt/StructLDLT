function SizeRankVersusSpeedup()
    # table of n = 500, 1000, 1500, 2000; r = n, n/2, n/10; f = 0, 1
    # ratio of [L,D,P] = UnStructCentro(A) divided by [Lp,Dp,Pp,Lm,Dm,Pm] = StructCentro(A)

    tol = 1e-6

    n_range = [500, 1000, 1500, 2000]
#    n_range = [1500, 1500, 1500, 1500]
    unstruct_times = zeros(4,3,2)
    struct_times = zeros(4,3,2)

    for n = 1:length(n_range)
        r_range = [n, n/2, n/10]
        for r = 1:length(r_range)
            for f = [0,1]
                A = RandCentro(n_range[n],r_range[r],f)
                tic()
                L = UnStructCentro(A,tol,192)
                unstruct_times[n,r,f+1] = toc()
                tic()
                Lp,Lm = StructCentro(A,tol)
                struct_times[n,r,f+1] = toc()
            end
        end
    end

display(unstruct_times)
display(struct_times)
display(unstruct_times./struct_times)

#    display(unstruct_times)
#    println()
#    display(mean(
#                 filter(
#                     x -> x < 1 ? 0 : x,
#                     unstruct_times),[1 2]
#                 )
#            )
#    println()
#    display(struct_times)
#    println()
#    display(mean(
#                 filter(
#                     x -> x < 1 ? 0.0 : x,struct_times),[1 2]
#                 )
#            )
#    println()
#    display(unstruct_times./struct_times)
#    println()
#    println(mean(
#                 filter(
#                     x -> x > 1,vec(unstruct_times[:,:,1])./vec(struct_times[:,:,1]))
#                 )
#            )
#    println(mean(
#                 filter(
#                    x -> x > 1,vec(unstruct_times[:,:,2])./vec(struct_times[:,:,2]))
#                 )
#            )
end

function RandCentro(n,r,f)
    # Returns a rank-r matrix with Centrosymmetry
    A = zeros(n,n)
    for i=1:r
        v = randn(n,1)
        if f == 1 && mod(i,2) == 0
            v -= v[n:-1:1]
        else
            v += v[n:-1:1]
        end
        A += v*v'
    end
    return A
end

function UnStructCentro(A,tol,block_size)
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
        G = lapack_full_fact(D[piv[j]],tol)
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


#function UnStructCentro(A,tol)
#    # Computes rank-revealing Cholesky factorization of A
#    Aout, piv, rank, info = LAPACK.pstrf!('L', A, tol)
##    L = tril(Aout)
#    L[piv,1:rank] = L[:,1:rank]
#    return L[:,1:rank]
#    return Aout
#end

function StructCentro(A,tol)
    # Computes the rank-revealing Cholesky factorization of A, utilizing Centrosymmetric structure
    n = size(A,1)
#    if mod(n,2) == 0
    m = convert(Int64,n/2)
    A11 = A[1:m,1:m]
    A12 = A[1:m,m+1:n]
#        Aout = 0
#        E = eye(m,m)[:,m:-1:1]
#    B11 = A11 + A12[:,m:-1:1]
#    Aout, piv, rank1, info = LAPACK.pstrf!('L', A11, tol)
    L1 = UnStructCentro(A11,tol,192)
#        Aout, piv, rank1, info = LAPACK.pstrf!('L', A11, tol)
#        Aout, piv, rank1, info = LAPACK.pstrf!('L', A12, tol)
        
#        L1 = tril(Aout)
#        L1[piv,1:rank1] = L1[:,1:rank1]
#    B22 = A11 - A12[:,m:-1:1]
    L2 = UnStructCentro(A12,tol,192)
#    Aout, piv, rank2, info = LAPACK.pstrf!('L', A12, tol)
#        L2 = tril(Aout)
#        L2[piv,1:rank2] = L2[:,1:rank2]
#    end
#    return L1[:,1:rank1], L2[:,1:rank2]
#    return Aout,Aout
    return 0,1
end

function TestRank()
    ranks = zeros(4,3,2)
    n_range = [500, 1000, 1500, 2000]
    for n = 1:4
        r_range = [n, n/2, n/10]
        for r = 1:3
            A = RandCentro(n_range[n],r_range[r])
            m = convert(Int64,n_range[n]/2)
            A11 = A[1:m,1:m]
            A12 = A[1:m,m+1:n_range[n]]
            B11 = A11 + A12[:,m:-1:1]
            display(B11)
            ranks[n,r,1] = rank(B11)
            B22 = A11 - A12[:,m:-1:1]
            display(B22)
            ranks[n,r,2] = rank(B22)
        end
    end
    display(ranks)
end

SizeRankVersusSpeedup()