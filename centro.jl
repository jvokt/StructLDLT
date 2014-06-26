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
                L = UnStructCentro(A,tol)
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

function UnStructCentro(A,tol)
    # Computes rank-revealing Cholesky factorization of A
    Aout, piv, rank, info = LAPACK.pstrf!('L', A, tol)
#    L = tril(Aout)
#    L[piv,1:rank] = L[:,1:rank]
#    return L[:,1:rank]
    return Aout
end

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
    Aout, piv, rank1, info = LAPACK.pstrf!('L', A11, tol)
#        Aout, piv, rank1, info = LAPACK.pstrf!('L', A11, tol)
#        Aout, piv, rank1, info = LAPACK.pstrf!('L', A12, tol)
        
#        L1 = tril(Aout)
#        L1[piv,1:rank1] = L1[:,1:rank1]                
#    B22 = A11 - A12[:,m:-1:1]
    Aout, piv, rank2, info = LAPACK.pstrf!('L', A12, tol)
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