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
        m = convert(Int64,n/2)
        r_range = [(m,m), (m,n/10), (m,0), (n/100,n/100), (n/100,0)]
        for r_i = 1:length(r_range)
            r1,r2 = convert((Int64,Int64),r_range[r_i])
            r = r1+r2
            println("n: ",n,", r1: ",r1,", r2: ",r2)
            B,C = RandCentro(n,r1,r2)
            A = FullCentro(B,C)
            for t_i = 1:trials
                Atemp = deepcopy(A)
                Btemp = deepcopy(B)
                Ctemp = deepcopy(C)

                tic() # <- start
                L,piv,rank,_ = UnStructCentro(Atemp)
                unstruct_times[n_i,r_i,t_i] = toc() # <- stop
                # L[piv,:] = tril(L)
                # L = L[:,1:rank]
                # println("error: ",norm(A-L*L'))
                
                tic() # <- start
                Lsym,psym,rsym,Lskew,pskew,rskew = StructCentro(Btemp,Ctemp)
                struct_times[n_i,r_i,t_i] = toc() # <- stop
                # Lsym[psym,:] = tril(Lsym)
                # Gsym = Lsym[:,1:rsym]
                # Gsym = [Gsym; Gsym[end:-1:1,:]]/sqrt(2)
                # Lskew[pskew,:] = tril(Lskew)
                # Gskew = Lskew[:,1:rskew]
                # Gskew = [Gskew; -Gskew[end:-1:1,:]]/sqrt(2)
                # println("error: ",norm(A-Gsym*Gsym'-Gskew*Gskew'))
            end
        end
    end
    display(mean(unstruct_times,3)./mean(struct_times,3))
    println()
    display(mean(unstruct_times./struct_times,3))
end

function CentroTrial1(n,m,r1,r2)
    A = randn(n,r)
    A = A*A'
    B1 = randn(m,r1)
    B1 = B1*B1'
    B2 = randn(m,r2)
    B2 = B2*B2'
    tic() # <- start
    # L = UnStructCentro(Acopy,tol)
    LAPACK.pstrf!('L', A, tol)
    unstruct_times[n_i,r_i,t_i] = toc() # <- stop
    # println("error: ",norm(A-L*L'))
    
    tic() # <- start
    # L1,L2 = StructCentro(A,tol)
    LAPACK.pstrf!('L', B1, tol)
    LAPACK.pstrf!('L', B2, tol)
    struct_times[n_i,r_i,t_i] = toc() # <- stop
    # println("error: ",norm(A-G1*G1'-G2*G2'))
end

function RandCentro(n,r1,r2)
    # Returns a rank-r matrix with Centrosymmetry
    m = convert(Int64,n/2)
    B = randn(m,r1)
    B = B*B'
    C = randn(m,r2)
    C = C*C'
    return B,C
end

function lapack_chol(A,tol=1e-6)
    A, piv, rank, info = LAPACK.pstrf!('L', A, tol)
    A[piv,:] = tril(A)
    return A[:,1:rank]
end

function FullCentro(B,C)
    X = (1/2)*(B + C)
    Y = (1/2)*(B - C)
    return [X Y[:,end:-1:1]; Y[end:-1:1,:] X[end:-1:1,end:-1:1]]
end

function UnStructCentro(A,tol=1e-6)
    # Computes the rank-revealing Cholesky factorization of A, without utilizing Centrosymmetric 
    return LAPACK.pstrf!('L', A, tol)
end

function UnStructCentro2(A,B,C,tol=1e-6)
    # Computes the rank-revealing Cholesky factorization of A, without utilizing Centrosymmetric structure
    #    return LAPACK.pstrf!('L', A, tol)
    X = (1/2)*(B + C)
    Y = (1/2)*(B - C)
    G11 = lapack_chol(X)
    X = (1/2)*(B + C)
    G21 = Y[end:-1:1,:]/G11'
    A22 = X[end:-1:1,end:-1:1] - G21*G21'
    G22 = lapack_chol(A22)
    G = [G11 zeros(size(G22)); G21 G22]
    return G
end

function StructCentro(B,C,tol=1e-6)
    # Computes the rank-revealing Cholesky factorization of A, utilizing Centrosymmetric structure
    Asym,psym,rsym,_ = LAPACK.pstrf!('L', B, tol)
    Askew,pskew,rskew,_ = LAPACK.pstrf!('L', C, tol)
    return Asym,psym,rsym,Askew,pskew,rskew
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