function SizeRankVersusSpeedup()
    # produce a table of n = 500, 1000, 1500, 2000
    # versus (r1,r2) = (m,m), (m,n/10), (m,0), (n/100,n/100), (n/100,0)
    # ratio of L,piv,rank = UnStructCentro(A) time divided
    # by Lsym,psym,rsym,Lskew,pskew,rskew = StructCentro(A) time

    tol = 1e-6
    trials = 5#0
    n_range = [1500, 3000, 4500, 6000]
    unstruct_times = zeros(length(n_range),2,trials)
    setup_times = zeros(length(n_range),2,trials)
    fact_times = zeros(length(n_range),2,trials)
    for n_i = 1:length(n_range)
        n = n_range[n_i]
        m = convert(Int64,n/2)
        r_range = [(m,m), (n/100,n/100)]
        for r_i = 1:length(r_range)
            r1,r2 = convert((Int64,Int64),r_range[r_i])
            r = r1+r2
            println("n: ",n,", r1: ",r1,", r2: ",r2)
            B,C = RandCentro(n,r1,r2)
            A = FullCentro(B,C)
            for t_i = 1:trials
                Atemp = deepcopy(A)
                tic() # <- start
                L,piv,rank,_ = UnStructCentro(Atemp)
                unstruct_times[n_i,r_i,t_i] = toc() # <- stop
                tic() # <- start
                X = A[1:m,1:m]
                Y = A[1:m,end:-1:m+1]
                B = X + Y
                C = X - Y
                setup_times[n_i,r_i,t_i] = toc()
                tic()
                Asym,psym,rsym,_ = LAPACK.pstrf!('L', B, tol)
                Askew,pskew,rskew,_ = LAPACK.pstrf!('L', C, tol)
                fact_times[n_i,r_i,t_i] = toc() # <- stop
            end
        end
    end
    struct_times = setup_times + fact_times
    #    display(mean(unstruct_times,3)./mean(struct_times,3))
    #    println()
    #    display(mean(unstruct_times./struct_times,3))
    #    println()
    #    display(mean(fact_times,3)./mean(struct_times,3))
    for n_i=1:length(n_range)
        Tu1 = mean(unstruct_times[n_i,1,:])
        Ts1 = mean(struct_times[n_i,1,:])
        Tsetup1 = mean(setup_times[n_i,1,:])
        Tu2 = mean(unstruct_times[n_i,2,:])
        Ts2 = mean(struct_times[n_i,2,:])
        Tsetup2 = mean(setup_times[n_i,2,:])
        t = [Tu1/Ts1, Tsetup1/Ts1, Tu2/Ts2, Tsetup2/Ts2]
        println(n_range[n_i]," & ",t[1]," & ",t[2],
                " & ",t[3]," & ",t[4],"  \\rule[-4pt]{0pt}{14pt}\\\\")
    end
#    n_range = [1500, 3000, 4500, 6000]
    speedup1 = zeros(length(n_range),2)
    for n_i = 1:length(n_range)
        n = n_range[n_i]
        m = convert(Int64,n/2)
        r_range = [(m,m), (n/100,n/100)]
        for r_i = 1:length(r_range)
            r1,r2 = convert((Int64,Int64),r_range[r_i])
            r = r1+r2
#            println("n: ",n,", r1: ",r1,", r2: ",r2)
            speedup1[n_i,r_i] = 2r^2/(2n+r1^2+r2^2)
        end
    end
    #display(speedup)
    #println()
    speedup2 = zeros(length(n_range),2)
    for n_i = 1:length(n_range)
        n = n_range[n_i]
        m = convert(Int64,n/2)
        r_range = [(m,m), (n/100,n/100)]
        for r_i = 1:length(r_range)
            r1,r2 = convert((Int64,Int64),r_range[r_i])
            r = r1+r2
#            println("n: ",n,", r1: ",r1,", r2: ",r2)
            speedup2[n_i,r_i] = 2n/(2n+r1^2+r2^2)
        end
    end
    #display(speedup)
    #println()
    for n_i=1:length(n_range)
#        t = [Tu1/Ts1, Tsetup1/Ts1, Tu2/Ts2, Tsetup2/Ts2]
        t = [speedup1[n_i,1], speedup2[n_i,1], 
             speedup1[n_i,2], speedup2[n_i,2]]
        println(n_range[n_i]," & ",t[1]," & ",t[2],
                " & ",t[3]," & ",t[4],"  \\rule[-4pt]{0pt}{14pt}\\\\")
    end
end

function UnStructCheck(L,piv,rank)
    L[piv,:] = tril(L)
    L = L[:,1:rank]
    println("error: ",norm(A-L*L'))
end

function StructCheck(Lsym,psym,rsym,Lskew,pskew,rskew)
    Lsym[psym,:] = tril(Lsym)
    Gsym = Lsym[:,1:rsym]
    Gsym = [Gsym; Gsym[end:-1:1,:]]/sqrt(2)
    Lskew[pskew,:] = tril(Lskew)
    Gskew = Lskew[:,1:rskew]
    Gskew = [Gskew; -Gskew[end:-1:1,:]]/sqrt(2)
    println("error: ",norm(A-Gsym*Gsym'-Gskew*Gskew'))    
end

function QE(m)
    I = speye(m)/sqrt(2)
    
    Qsym = spzeros(2*m,m)
    Qsym[1:m,1:m] = I
    Qsym[m+1:end,1:m] = I[m:-1:1,:]
    
    Qskew = spzeros(2*m,m)
    Qskew[1:m,1:m] = I
    Qskew[m+1:end,1:m] = -I[m:-1:1,:]
    return Qsym,Qskew
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


function CompareStructCentro()
    n_range = [500 1000 1500 2000]
    unstruct_times = zeros(length(n_range),5)
    struct_times1 = zeros(length(n_range),5)
    struct_times2 = zeros(length(n_range),5)
    struct_times3 = zeros(length(n_range),5)
    struct_times4 = zeros(length(n_range),5)
    for n_i = 1:length(n_range)
        n = n_range[n_i]
        m = int(n/2)
        Qsym,Qskew = QE(m)
        r_range = [(m,m), (m,n/10), (m,0), (n/100,n/100), (n/100,0)]
        for r_i = 1:length(r_range)
            r1,r2 = convert((Int64,Int64),r_range[r_i])
            r = r1+r2
            println("n: ",n,", r1: ",r1,", r2: ",r2)
            B,C = RandCentro(n,r1,r2)
            A = FullCentro(B,C)
            Atemp = deepcopy(A)
            
            tic()
            UnStructCentro(Atemp)
            unstruct_times[n_i,r_i] = toc()
            
            tic()
            StructCentro1(A,Qsym,Qskew)
            struct_times1[n_i,r_i] = toc()

            tic()
            StructCentro2(A)
            struct_times2[n_i,r_i] = toc()
            
            X = A[1:end/2,1:end/2]
            Y = A[1:end/2,end:-1:end/2+1]
            tic()
            StructCentro3(X,Y)
            struct_times3[n_i,r_i] = toc()
            
            tic()
            StructCentro4(B,C)
            struct_times4[n_i,r_i] = toc()
        end
    end
    display(unstruct_times./struct_times1)
    println()
    display(unstruct_times./struct_times2)
    println()
    display(unstruct_times./struct_times3)
    println()
    display(unstruct_times./struct_times4)
    println()
end

function StructCentro1(A,Qsym,Qskew,tol=1e-6)
    # Computes the rank-revealing Cholesky factorization of A, utilizing Centrosymmetric structure
    B = Qsym'*A*Qsym
    Asym,psym,rsym,_ = LAPACK.pstrf!('L', B, tol)
    C = Qskew'*A*Qskew
    Askew,pskew,rskew,_ = LAPACK.pstrf!('L', C, tol)
    return Asym,psym,rsym,Askew,pskew,rskew
end

function StructCentro2(A,tol=1e-6)
    # Computes the rank-revealing Cholesky factorization of A, utilizing Centrosymmetric structure
    X = A[1:end/2,1:end/2]
    Y = A[1:end/2,end:-1:end/2+1]
    B = X + Y
    C = X - Y
    Asym,psym,rsym,_ = LAPACK.pstrf!('L', B, tol)
    Askew,pskew,rskew,_ = LAPACK.pstrf!('L', C, tol)
    return Asym,psym,rsym,Askew,pskew,rskew
end

function StructCentro3(X,Y,tol=1e-6)
    # Computes the rank-revealing Cholesky factorization of A, utilizing Centrosymmetric structure
#    X = A[1:end/2,1:end/2]
#    Y = A[1:end/2,end:-1:end/2+1]
    B = X + Y
    Asym,psym,rsym,_ = LAPACK.pstrf!('L', B, tol)
    C = X - Y
    Askew,pskew,rskew,_ = LAPACK.pstrf!('L', C, tol)
    return Asym,psym,rsym,Askew,pskew,rskew
end

function StructCentro4(B,C,tol=1e-6)
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

function EstimateConstant()
    n_range = [500, 1000, 1500, 2000, 3000]
    times = zeros(length(n_range))
    for n_i = 1:length(n_range)
        n = n_range[n_i]
        A = randn(n,n); A = A*A'
        tic()
        LAPACK.pstrf!('L',A,1e-6)
        times[n_i] = toc()
    end
    figure()
    plot(n_range,times)
    figure()
    plot(log(n_range),log(times))
    c = [log(n_range) ones(length(n_range))]\log(times)
end

function PredictedSpeedup()
    c = 1/2 # EstimateConstant()
#    n_range = [500, 1000, 1500, 2000]
    n_range = [1500, 3000, 4500, 6000]
    speedup = zeros(length(n_range),5)
    for n_i = 1:length(n_range)
        n = n_range[n_i]
        m = convert(Int64,n/2)
        r_range = [(m,m), (m,n/10), (m,0), (n/100,n/100), (n/100,0)]
        for r_i = 1:length(r_range)
            r1,r2 = convert((Int64,Int64),r_range[r_i])
            r = r1+r2
            println("n: ",n,", r1: ",r1,", r2: ",r2)
            speedup[n_i,r_i] = 2r^2/(2n+r1^2+r2^2)
        end
    end
    display(speedup)
    println()
end

SizeRankVersusSpeedup()
#TestRank()
#CompareStructCentro()
#PredictedSpeedup()