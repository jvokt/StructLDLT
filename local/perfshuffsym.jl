function SizeRankVersusSpeedup()
    # produce a table of n = 500, 1000, 1500, 2000
    # versus (r1,r2) = (m,m), (m,n/10), (m,0), (n/100,n/100), (n/100,0)
    # ratio of L,piv,rank = UnStructPerfShuff(A) time divided
    # by Lsym,psym,rsym,Lskew,pskew,rskew = StructPerfShuff(A) time    
    tol = 1e-6
    trials = 10
    num_ranks = 2
    n_range = [39, 55, 67, 77]
    unstruct_times = zeros(length(n_range),num_ranks)
    setup_times = zeros(length(n_range),num_ranks)
    fact_times = zeros(length(n_range),num_ranks)
    for n_i = 1:length(n_range)
        n = n_range[n_i]
        nsym = convert(Int64,n*(n+1)/2)
        nskew = convert(Int64,n*(n-1)/2)
        r_range = [(nsym,nskew), (n,n)]
        for r_i = 1:length(r_range)
            r1,r2 = convert((Int64,Int64),r_range[r_i])
            r = r1+r2
            println("n: ",n,", r1: ",r1,", r2: ",r2)
            B,C = RandPerfShuff(n,r1,r2)
            A = FullPerfShuff(B,C,n)
            As = [deepcopy(A) for i=1:trials]
            Bs = [deepcopy(B) for i=1:trials]
            Cs = [deepcopy(C) for i=1:trials]
            piv = zeros(n^2)
            rank = 0
            tic()
            for t_i = 1:trials
                L,piv,rank,_ = LAPACK.pstrf!('L', As[t_i], tol)
            end
            unstruct_times[n_i,r_i] = toc()
            UnStructCheck(A,As[1],piv,rank)

            tic()
            for t_i = 1:trials
                PS = PerfShuff(n,n)
                sym_indices = SymIndices(n)
                X = A[sym_indices,sym_indices]
                Y = A[sym_indices,PS[sym_indices]]
                B = X + Y
                skew_indices = SkewIndices(n)
                X = A[skew_indices,skew_indices]
                Y = A[skew_indices,PS[skew_indices]]
                C = X - Y
            end
            setup_times[n_i,r_i] = toc()

            psym = zeros(n^2)
            rsym = 0
            pskew = zeros(n^2)
            rskew = 0
            tic()
            for t_i = 1:trials
                Lsym,psym,rsym,_ = LAPACK.pstrf!('L', Bs[t_i], tol)
                Lskew,pskew,rskew,_ = LAPACK.pstrf!('L', Cs[t_i], tol)
            end
            fact_times[n_i,r_i] = toc()
            StructCheck(A,Bs[1],psym,rsym,Cs[1],pskew,rskew)
        end
    end
    struct_times = setup_times + fact_times
    for n_i=1:length(n_range)
        Tu1 = mean(unstruct_times[n_i,1,:])
        Ts1 = mean(struct_times[n_i,1,:])
        Tsetup1 = mean(setup_times[n_i,1,:])
        Tu2 = mean(unstruct_times[n_i,2,:])
        Ts2 = mean(struct_times[n_i,2,:])
        Tsetup2 = mean(setup_times[n_i,2,:])
        t = round([Tu1/Ts1, Tsetup1/Ts1, Tu2/Ts2, Tsetup2/Ts2],2)
        println(n_range[n_i]," & ",t[1]," & ",t[2],
                " & ",t[3]," & ",t[4],"  \\rule[-4pt]{0pt}{14pt}\\\\")
    end
end

function DevectorizedSetup()
#                B = zeros(nsym,nsym)
#                C = zeros(nskew,nskew)
                #  sympacked(i,j) = i-(1/2)(j-1)j+(j-1)n
                #  skewpacked(i,j) = i-(1/2)(j-1)j-j+(j-1)n
#                for j=1:n
#                    for i=j:n
#                        for l=1:n
#                            for k=l:n
#                                p = i+(j-1)n
#                                q = k+(l-1)n
#                                rsym = i-(1/2)*(j-1)j+(j-1)n
#                                rskew = i-(1/2)*(j-1)j-j+(j-1)n
#                                ssym = k-(1/2)*(l-1)l+(l-1)n
#                                sskew = k-(1/2)*(l-1)l-l+(l-1)n
#                                if i > j && k > l
#                                    # C
#                                    C[rskew,sskew] = A[p,q]-A[p,PS[q]]
#                                end
                                # B
#                                B[rsym,ssym] = A[p,q]+A[p,PS[q]]
#                            end
#                        end
#                    end
#                end
end

function PredictedSpeedups()
    n_range = [39, 55, 67, 77]
    num_ranks = 2
    speedup1 = zeros(length(n_range),num_ranks)
    for n_i = 1:length(n_range)
        n = n_range[n_i]
        nsym = convert(Int64,n*(n+1)/2)
        nskew = convert(Int64,n*(n-1)/2)
        r_range = [(nsym,nskew), (n,n)]
        for r_i = 1:length(r_range)
            r1,r2 = convert((Int64,Int64),r_range[r_i])
            r = r1+r2
#            println("n: ",n,", r1: ",r1,", r2: ",r2)
            speedup1[n_i,r_i] = 2n*r^2/(2n*(n^2+1)+(n+1)*r1^2+(n-1)r2^2)
        end
    end
    display(speedup1)
    println()
    speedup2 = zeros(length(n_range),num_ranks)
    for n_i = 1:length(n_range)
        n = n_range[n_i]
        nsym = convert(Int64,n*(n+1)/2)
        nskew = convert(Int64,n*(n-1)/2)
        r_range = [(nsym,nskew), (n,n)]
        for r_i = 1:length(r_range)
            r1,r2 = convert((Int64,Int64),r_range[r_i])
            r = r1+r2
#            println("n: ",n,", r1: ",r1,", r2: ",r2)
            speedup2[n_i,r_i] = 2n*(n^2+1)/(2n*(n^2+1)+(n+1)r1^2+(n-1)r2^2)
        end
    end
    display(speedup2)
    println()
    speedup3 = zeros(length(n_range),num_ranks)
    for n_i = 1:length(n_range)
        n = n_range[n_i]
        nsym = convert(Int64,n*(n+1)/2)
        nskew = convert(Int64,n*(n-1)/2)
        r_range = [(nsym,nskew), (n,n)]
        for r_i = 1:length(r_range)
            r1,r2 = convert((Int64,Int64),r_range[r_i])
            r = r1+r2
#            println("n: ",n,", r1: ",r1,", r2: ",r2)
            speedup3[n_i,r_i] = 2n*r^2/((n+1)*r1*(r1+2)+(n-1)*r2*(r2+2))
        end
    end
    display(speedup3)
    println()
    speedup4 = zeros(length(n_range),num_ranks)
    for n_i = 1:length(n_range)
        n = n_range[n_i]
        nsym = convert(Int64,n*(n+1)/2)
        nskew = convert(Int64,n*(n-1)/2)
        r_range = [(nsym,nskew), (n,n)]
        for r_i = 1:length(r_range)
            r1,r2 = convert((Int64,Int64),r_range[r_i])
            r = r1+r2
#            println("n: ",n,", r1: ",r1,", r2: ",r2)
            speedup4[n_i,r_i] = 2((n+1)r1+(n-1)r2)/((n+1)*r1*(r1+2)+(n-1)*r2*(r2+2))
        end
    end
    display(speedup4)
    println()
    for n_i=1:length(n_range)
        t = round([speedup1[n_i,1], speedup2[n_i,1],
             speedup1[n_i,2], speedup2[n_i,2]],2)
        println(n_range[n_i]," & ",t[1]," & ",t[2],
                " & ",t[3]," & ",t[4],"  \\rule[-4pt]{0pt}{14pt}\\\\")
    end
    for n_i=1:length(n_range)
        t = round([speedup3[n_i,1], speedup4[n_i,1],
             speedup3[n_i,2], speedup4[n_i,2]],2)
        println(n_range[n_i]," & ",t[1]," & ",t[2],
                " & ",t[3]," & ",t[4],"  \\rule[-4pt]{0pt}{14pt}\\\\")
    end
end

function UnStructCheck(A,L,piv,rank)
    L[piv,:] = tril(L)
    L = L[:,1:rank]
    println("error: ",norm(A-L*L'))
end

function StructCheck(A,Lsym,psym,rsym,Lskew,pskew,rskew)
    # n = 3;r1 = 6; r2 = 3
    # A = RandPerfShuff(n,r1,r2);
    # Lsym,psym,rsym,Lskew,pskew,rskew = StructPerfShuff3(A);
    # StructCheck(A,Lsym,psym,rsym,Lskew,pskew,rskew)
    n2 = size(A,1)
    n = int(sqrt(n2))
#    nsym = int(n*(n+1)/2)
#    nskew = int(n*(n-1)/2)
#    s = snn(n)
#    Qsym,Qskew = Qnn(n)
    Lsym[psym,:] = tril(Lsym)
    #Gsym = scale(s,Lsym[:,1:rsym])
    #Gsym = Qsym*Gsym
    Gsym = Lsym[:,1:rsym]
    
#    Gsym = SymUnpackRows2(Lsym[:,1:rsym],n)
    Lskew[pskew,:] = tril(Lskew)
    Gskew = Lskew[:,1:rskew]
#    Gskew = Qskew*Gskew
    Ar = FullPerfShuff(Gsym*Gsym',Gskew*Gskew',n)
    println("error: ",norm(A-Ar))
end

function Qnn(n)
    # Q is the n^2 by n^2 orthgonal matrix which 
    # block diagonalizes symmetric perfect shuffle invariant matrices

    # Define Qsym
    Qsym = spzeros(n^2,int(n*(n+1)/2))
    e = eye(n,n)
    k = 1
    for i=1:n
        for j=i:n
            if i == j
                Qsym[i+(i-1)*n,k] = 1
            else
                Qsym[i+(j-1)*n,k] = 1/sqrt(2)
                Qsym[j+(i-1)*n,k] = 1/sqrt(2)
            end
            k = k+1
        end
    end
    # Define Qskew
    Qskew = spzeros(n^2,int(n*(n-1)/2))
    k = 1
    for i=1:n
        for j=i+1:n
            Qskew[j+(i-1)*n,k] = 1/sqrt(2)
            Qskew[i+(j-1)*n,k] = -1/sqrt(2)
            k = k+1
        end
    end
    return Qsym,Qskew
end

function Onn(n)
    # O is the n^2 by n^2 zero-one matrix which 
    # block diagonalizes symmetric perfect shuffle invariant matrices

    # Define Osym
    Osym = spzeros(n^2,int(n*(n+1)/2))
    e = eye(n,n)
    k = 1
    for i=1:n
        for j=i:n
            if i == j
                Osym[i+(i-1)*n,k] = 1
            else
                Osym[i+(j-1)*n,k] = 1
                Osym[j+(i-1)*n,k] = 1
            end
            k = k+1
        end
    end
    # Define Oskew
    Oskew = spzeros(n^2,int(n*(n-1)/2))
    k = 1
    for i=1:n
        for j=i+1:n
            Oskew[j+(i-1)*n,k] = 1
            Oskew[i+(j-1)*n,k] = -1
            k = k+1
        end
    end
end

function tnn(n)
    unity = [j+(j-1)*n-(j-1)j/2 for j=1:n]
    t = ones(n^2)/sqrt(2)
    t[unity] = 1
    return t
end

function snn(n)
    unity = [j+(j-1)*n-(j-1)j/2 for j=1:n]
    t = ones(int(n*(n+1)/2))
    t[unity] = 1/sqrt(2)
    return t
end

function SymUnpackRows1(B,n)
    # [unpackIndex(c,n) for c=1:n*(n+1)/2]
    G = zeros(n^2,size(B,2))
    for c=1:n*(n+1)/2
        i, j = unpackIndex(c,n)
        if i > j
            r1 = i+(j-1)*n
            G[r1,:] = B[c,:]
            r2 = j+(i-1)*n
            G[r2,:] = B[c,:]
        else
            r1 = i+(j-1)*n
            G[r1,:] = B[c,:]
        end            
    end
    return G
end

function SymUnpackRows2(B,n)
    G = zeros(n^2,size(B,2))
    sym_indices = SymIndices(n)
    PS = PerfShuff(n,n)
    G[sym_indices,:] = B
    G[PS[sym_indices],:] = B
    return G
end

function SymIndices(n)
    nsym = int(n*(n+1)/2)
    sym_indices = zeros(Int64,nsym)
    idx = 1
    for j=1:n
        for i=j:n
            sym_indices[idx] = i+(j-1)*n
            idx += 1
        end
    end
    return sym_indices
end

function SkewIndices(n)
    nskew = int(n*(n-1)/2)
    skew_indices = zeros(Int64,nskew)
    idx = 1
    for j=1:n
        for i=j+1:n
            skew_indices[idx] = i+(j-1)*n
            idx += 1
        end
    end
    return skew_indices
end

function UnpackSkewIndex(c,n)
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

function SkewUnpackRows(B,n)
    # [unpackIndex(c,n) for c=1:n*(n+1)/2]
    G = zeros(n^2,size(B,2))
    for c=1:n*(n-1)/2
        i, j = skewIndex(c,n)
        if i > j
            r1 = i+(j-1)*n
            G[r1,:] = B[c,:]#/sqrt(2)
            r2 = j+(i-1)*n
            G[r2,:] = B[c,:]#/sqrt(2)
        else
            r1 = i+(j-1)*n
            G[r1,:] = B[c,:]
        end            
    end
    return G
end

function CompareSymUnpackRows()
    n_range = [20, 40, 60, 80, 100]
    for n_i = 1:length(n_range)
        n = n_range[n_i]
        nsym = convert(Int64,n*(n+1)/2)
        B = randn(nsym,6n)
        Qsym,Qskew = Qnn(n)
        Osym,Oskew = Onn(n)
        
        tic()
        G = Osym*B
        toc()
        
        tic()
        G1 = SymUnpackRows1(B,n)
        toc()
        
        tic()
        G2 = SymUnpackRows2(B,n)
        toc()
        
        println("Error: ",norm(G-G1))
        println("Error: ",norm(G-G2))
    end
end

function CompareSymBlock()
    n_range = [20, 40, 60]#, 80, 100]
    for n_i = 1:length(n_range)
        n = n_range[n_i]
        nsym = convert(Int64,n*(n+1)/2)
        nskew = convert(Int64,n*(n-1)/2)
        Qsym,Qskew = Qnn(n)
        r_range = [(nsym,nskew), (nsym,n/2), (nsym,0), (n/2,n/2), (n,0)]
        for r_i = 1:length(r_range)
            r1,r2 = convert((Int64,Int64),r_range[r_i])
            println("n: ",n,", r1: ",r1,", r2: ",r2)
            A = RandPerfShuff(n,r1,r2)
            
            tic()
            B1 = Qsym'*A*Qsym
            toc()
            
            tic()
            B2 = sym_block(A)
            toc()

            tic()
            B3 = struct_block(A,1,1,nsym,nsym)
            toc()
            println("error: ",norm(B1-B2))
        end
    end
end

function CompareRand()
    n_range = [20, 40, 60]#, 80, 100]
    for n_i = 1:length(n_range)
        n = n_range[n_i]
        nsym = convert(Int64,n*(n+1)/2)
        nskew = convert(Int64,n*(n-1)/2)
        r_range = [(nsym,nskew), (nsym,n/2), (nsym,0), (n/2,n/2), (n,0)]
        for r_i = 1:length(r_range)
            r1,r2 = convert((Int64,Int64),r_range[r_i])
            println("n: ",n,", r1: ",r1,", r2: ",r2)
            
            tic()
            RandPerfShuff1(n,r1,r2)
            toc()
            
            tic()
            RandPerfShuff(n,r1,r2)
            toc()
        end
    end
end

function SkewInSymIndices(n)
    nskew = int(n*(n-1)/2)
    indices = zeros(Int64,nskew)
    idx = 1
    idx2 = 1
    for j=1:n
        for i=j:n
            if i > j
                indices[idx] = idx2
                idx += 1
            end
            idx2 += 1
        end
    end
    return indices
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

function FullPerfShuff(B,C,n)
    PS = PerfShuff(n,n)
    sym = SymIndices(n)
    skew = SkewIndices(n)
    
    A1 = zeros(n^2,n^2)
    A1[sym,sym] = B
    A1[PS[sym],PS[sym]] = B
    A1[sym,PS[sym]] = B
    A1[PS[sym],sym] = B
    
    A2 = zeros(n^2,n^2)
    A2[skew,skew] = C
    A2[PS[skew],PS[skew]] = C
    A2[skew,PS[skew]] = -C
    A2[PS[skew],skew] = -C
    return A1+A2
end

function RandPerfShuff1(n,r1,r2)
    # Returns a rank-r matrix with Perfect Shuffle Symmetry
    nsym = convert(Int64,n*(n+1)/2)
    nskew = convert(Int64,n*(n-1)/2)
    B = randn(nsym,r1)
    B = B*B'
    C = randn(nskew,r2)
    C = C*C'
    Qsym,Qskew = Qnn(n)
    X = Qsym*B*Qsym'
    Y = Qskew*C*Qskew'
    return X + Y
end

function RandPerfShuff2(n,r1,r2)
    # Returns a rank-r matrix with Perfect Shuffle Symmetry
    nsym = convert(Int64,n*(n+1)/2)
    nskew = convert(Int64,n*(n-1)/2)
    Qsym,Qskew = Qnn(n)
    B = Qsym*randn(nsym,r1)
    X = B*B'
    C = Qskew*randn(nskew,r2)
    Y = C*C'
    return X + Y
end

function RandPerfShuff3(n,r1,r2)
    # Returns a rank-r matrix with Perfect Shuffle Symmetry
    nsym = convert(Int64,n*(n+1)/2)
    nskew = convert(Int64,n*(n-1)/2)
    B = SymUnpackRows(randn(nsym,r1),n)
    X = B*B'
    C = SkewUnpackRows(randn(nskew,r2),n)
    Y = C*C'
#    Qsym,Qskew = Qnn(n)
#    X = Qsym*B*Qsym'
#    Y = Qskew*C*Qskew'
    return X + Y
end

function lapack_chol(A,tol=1e-6)
    A, piv, rank, info = LAPACK.pstrf!('L', A, tol)
    A[piv,:] = tril(A)
    return A[:,1:rank]
end

function FullPerfShuff4(B,C)
    Qsym,Qskew = Qnn(n)
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

function StructPerfShuff1(A,tol=1e-6)
    # Computes the rank-revealing Cholesky factorization of A, utilizing PerfShuffSymmetric structure
    n2 = size(A,1)
    n = int(sqrt(n2))
    nsym = int(n*(n+1)/2)
    Qsym,Qskew = Qnn(n)
    B = Qsym'*A*Qsym
    C = Qskew'*A*Qskew
    Asym,psym,rsym,_ = LAPACK.pstrf!('L', B, tol)
    Askew,pskew,rskew,_ = LAPACK.pstrf!('L', C, tol)
    return Asym,psym,rsym,Askew,pskew,rskew
end

function StructPerfShuff2(A,tol=1e-6)
    # Computes the rank-revealing Cholesky factorization of A, utilizing PerfShuffSymmetric structure
    n2 = size(A,1)
    n = int(sqrt(n2))
#    tic()
    B = sym_block(A,n)
#    toc()
#    tic()
    C = skew_block(A,n)
#    toc()
#    tic()
    Asym,psym,rsym,_ = LAPACK.pstrf!('L', B, tol)
#    toc()
#    tic()
    Askew,pskew,rskew,_ = LAPACK.pstrf!('L', C, tol)
#    toc()
    return Asym,psym,rsym,Askew,pskew,rskew
end

function StructPerfShuff3(A,tol=1e-6)
    # Computes the rank-revealing Cholesky factorization of A, utilizing PerfShuffSymmetric structure
    n2 = size(A,1)
    n = int(sqrt(n2))
    PS = PerfShuff(n,n)
    nsym = int(n*(n+1)/2)
    nskew = int(n*(n-1)/2)
    sym_indices = zeros(Int64,nsym)
    skew_indices = zeros(Int64,nskew)
    sym_idx = 1
    skew_idx = 1
    for j=1:n
        j_base = (j-1)*n
        sym_indices[sym_idx] = j+j_base
        sym_idx += 1
        for i=j+1:n
            sym_indices[sym_idx] = i+j_base
            sym_idx += 1
            skew_indices[skew_idx] = i+j_base
            skew_idx += 1
        end
    end
    X = A[sym_indices,sym_indices]
    Y = A[sym_indices,PS[sym_indices]]
    # s = snn(n)
    # B = scale(s,scale(X + Y,s))
    B = X + Y
    X = A[skew_indices,skew_indices]
    Y = A[skew_indices,PS[skew_indices]]
    C = X - Y
    Asym,psym,rsym,_ = LAPACK.pstrf!('L', B, tol)
    Askew,pskew,rskew,_ = LAPACK.pstrf!('L', C, tol)
    return Asym,psym,rsym,Askew,pskew,rskew
end

function sym_block(A,n)
    PS = PerfShuff(n,n)
    sym_indices = SymIndices(n)
    X = A[sym_indices,sym_indices]
    Y = A[sym_indices,PS[sym_indices]]
    s = snn(n)
    return scale(s,scale(X+Y,s))
end

function skew_block(A,n)
    PS = PerfShuff(n,n)
    skew_indices = SkewIndices(n)
    X = A[skew_indices,skew_indices]
    Y = A[skew_indices,PS[skew_indices]]
    return X - Y
end

function sym_block_precompute(A)
    n2 = size(A,1)
    n = int(sqrt(n2))
    nsym = int(n*(n+1)/2)
    sym_indices = zeros(nsym)
    idx = 1
    for j=1:n
        for i=j:n
            sym_indices[idx] = i+(j-1)*n
        end
    end
    return A[sym_indices,sym_indices]
end

function StructPerfShuff4(B,C,tol=1e-6)
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
n = p*r
v = zeros(Int64,n)
for k=1:r
    v[(k-1)*p+1:k*p] = k:r:n
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

function CompareStructPerfShuff()
    n_range = [20, 40, 60]#, 80, 100]
    unstruct_times = zeros(length(n_range),5)
    struct_times1 = zeros(length(n_range),5)
    struct_times2 = zeros(length(n_range),5)
    struct_times3 = zeros(length(n_range),5)
#    struct_times4 = zeros(length(n_range),5)
    for n_i = 1:length(n_range)
        n = n_range[n_i]
        nsym = convert(Int64,n*(n+1)/2)
        nskew = convert(Int64,n*(n-1)/2)
        r_range = [(nsym,nskew), (nsym,n/2), (nsym,0), (n/2,n/2), (n,0)]
        for r_i = 1:length(r_range)
            r1,r2 = convert((Int64,Int64),r_range[r_i])
            r = r1+r2
            println("n: ",n,", r1: ",r1,", r2: ",r2)
            A = RandPerfShuff(n,r1,r2)
            Atemp = deepcopy(A)
            
            tic()
            UnStructPerfShuff(Atemp)
            unstruct_times[n_i,r_i] = toc()
            
            tic()
            StructPerfShuff1(A)
            struct_times1[n_i,r_i] = toc()

            tic()
            StructPerfShuff2(A)
            struct_times2[n_i,r_i] = toc()
            
            #X = A[1:end/2,1:end/2]
            #Y = A[1:end/2,end:-1:end/2+1]
            tic()
            StructPerfShuff3(A)
            struct_times3[n_i,r_i] = toc()
            
            #tic()
            #StructPerfShuff4(B,C)
            #struct_times4[n_i,r_i] = toc()
        end
    end
    display(unstruct_times./struct_times1)
    println()
    display(unstruct_times./struct_times2)
    println()
    display(unstruct_times./struct_times3)
    println()
    #display(unstruct_times./struct_times4)
    #println()
end

function PredictedSpeedup()
    n_range = [20, 40, 60, 80]#, 100]
    speedup = zeros(length(n_range),5)
    for n_i = 1:length(n_range)
        n = n_range[n_i]
        nsym = convert(Int64,n*(n+1)/2)
        nskew = convert(Int64,n*(n-1)/2)
        r_range = [(nsym,nskew), (nsym,6n), (nsym,0), (6n,6n), (6n,0)]
        for r_i = 1:length(r_range)
            r1,r2 = convert((Int64,Int64),r_range[r_i])
            r = r1+r2
            println("n: ",n,", r1: ",r1,", r2: ",r2)
            speedup[n_i,r_i] = 2n*r^2/(n*(n+1)^2+n*(n-1)^2+(n+1)r1^2+(n-1)r2^2)
        end
    end
    display(speedup)
    println()
end

SizeRankVersusSpeedup()
#TestRank()
#CompareSymUnpackRows()
#CompareRand()
#CompareSymBlock()
#CompareStructPerfShuff()
#PredictedSpeedups()