using PyPlot

function TestIntegralTransformation()

    molecules = ["HF","NH3","H2O2","N2H4","C2H5OH"]
    num_basis_fns = [34, 48, 68, 82, 123]
    nmols = 2#length(num_basis_fns)
    gen_plots = false
    tol = 1e-6
    
    unfact_time = zeros(nmols)
    unfact_space = zeros(nmols)
    unfact_error = zeros(nmols)
    unfact_fevals = zeros(nmols)

    full_fact_time = zeros(nmols)
    full_fact_space = zeros(nmols)
    full_fact_error = zeros(nmols)
    full_fact_fevals = zeros(nmols)

    full_trans_time = zeros(nmols)
    full_trans_space = zeros(nmols)
    full_trans_error = zeros(nmols)

    struct_fact_time = zeros(nmols)
    struct_fact_space = zeros(nmols)
    struct_fact_error = zeros(nmols)
    struct_fact_fevals = zeros(nmols)

    struct_trans_time = zeros(nmols)
    struct_trans_space = zeros(nmols)
    struct_trans_error = zeros(nmols)    
    
    for i=1:nmols
        molecule = molecules[i]
        n = num_basis_fns[i]
        r = ranks[i]
        A = RandERIMatrix(n,r)
        C = randn(n,n)
        
        println("Molecule: ",molecule)

        # Unfactorized transformation
        println("Unfactorized transformation")
        tic()
        M = unfact_trans(A,C)
        unfact_time[i] = toc()
        unfact_error[i] = norm(kron(C,C)*A*kron(C,C)'-M,1)
        println("error: ",unfact_error[i])

        # Full factorization
        println("Full factorization")
        tic()
        L, piv, full_fact_space[i], full_fact_fevals[i] = lapack_full_fact(A,tol)
        full_fact_time[i] = toc()
        full_fact_error[i] = norm(A[piv,piv] - L*L',1)
        println("error: ", full_fact_error[i])

        # Full transformation
        println("Full transformation")
        tic()
        B, full_trans_space[i] = full_trans(L,piv,C)
        full_trans_time[i] = toc()
        full_trans_error[i] = norm(kron(C,C)*A*kron(C,C)' - B*B',1)
        println("error: ",full_trans_error[i])

        println("Total time: ",full_fact_time[i]+full_trans_time[i])
        
        # Structured factorization
        println("Structured factorization")
        tic()
        L, piv, struct_fact_space[i], struct_fact_fevals[i] = lapack_struct_fact(A,tol)
        struct_fact_time[i] = toc()
        J = QsymTransform(L,piv,n)
        struct_fact_error[i] = norm(A - J*J',1)
        println("error: ",struct_fact_error[i])
        
        # Structured transformation
        println("Structured transformation")
        tic()
        B, struct_trans_space[i] = struct_trans(L,piv,C)
        struct_trans_time[i] = toc()
        struct_trans_error[i] = norm(kron(C,C)*A*kron(C,C)' - B*B',1)
        println("error: ",struct_trans_error[i])

        println("Total time: ",struct_fact_time[i]+struct_trans_time[i])
    end
    
    if gen_plots
        plotComparison(num_basis_fns[1:nmols],
                       full_fact_time,
                       "Full factorization",
                       struct_fact_time,
                       "Structured factorization",
                       "Number of basis functions",
                       "Time (seconds)","Time to factorize")
        plotComparison(num_basis_fns[1:nmols],
                       full_fact_space,
                       "Full factorization",
                       struct_fact_space,
                       "Structured factorization",
                       "Number of basis functions",
                       "Space (in Float64s)",
                       "Space to factorize")
        plotComparison(num_basis_fns[1:nmols],
                       full_fact_error,
                       "Full factorization",
                       struct_fact_error,
                       "Structured factorization",
                       "Number of basis functions",
                       "Error (1-norm)","Error to factorize")
        plotComparison(num_basis_fns[1:nmols],
                       full_fact_fevals,
                       "Full factorization",
                       struct_fact_fevals,
                       "Structured factorization",
                       "Number of basis functions",
                       "Fevals","Fevals to factorize")

        plotComparison(num_basis_fns[1:nmols],
                       full_trans_time,
                       "Full transformation",
                       struct_trans_time,
                       "Structured transformation",
                       "Number of basis functions",
                       "Time (seconds)",
                       "Time to transform")
        plotComparison(num_basis_fns[1:nmols],
                       full_trans_space,
                       "Full transformation",
                       struct_trans_space,
                       "Structured transformation",
                       "Number of basis functions",
                       "Space (in Float64s)","Space to transform")
        plotComparison(num_basis_fns[1:nmols],
                       full_trans_error,
                       "Full transformation",
                       struct_trans_error,
                       "Structured transformation",
                       "Number of basis functions",
                       "Error (1-norm)","Error to transform")

        plotComparison(num_basis_fns[1:nmols],
                       full_fact_time+full_trans_time,
                       "Full factorization plus transformation",
                       struct_fact_time+struct_trans_time,
                       "Structured factorization plus transformation",
                       "Number of basis functions",
                       "Time (seconds)",
                       "Time to factorize and transform")
        plotComparison(num_basis_fns[1:nmols],
                       full_fact_space+full_trans_space,
                       "Full factorization plus transformation",
                       struct_fact_space+struct_trans_space,
                       "Structured factorization plus transformation",
                       "Number of basis functions",
                       "Space in Float64s)",
                       "Space to factorize and transform")
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

function unfact_trans(A,C)
    # Returns the transformed integrals by a series of modal products
    n,n = size(C)
    B0 = deepcopy(A)
    B1 = C*reshape(A,n,n^3)
    B2 = C*reshape(permutedims(reshape(B1,n,n,n,n),[2 1 3 4]),n,n^3)
    B3 = C*reshape(permutedims(reshape(B2,n,n,n,n),[3 2 1 4]),n,n^3)
    B4 = C*reshape(permutedims(reshape(B3,n,n,n,n),[4 2 3 1]),n,n^3)
    B = reshape(permutedims(reshape(B4,n,n,n,n),[2 3 4 1]),n^2,n^2)
    return B
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

function blocked_full_fact(A,tol,block_size)
    # Computes the blocked full factorization of A
    full_fact_space = 0
    full_fact_fevals = 0
    n = size(A,1)
    num_blocks = convert(Int64,div(n,block_size)) + (convert(Bool,mod(n,block_size)) ? 1 : 0)
    L = Array{Float64,2}[]
    # evaluate diagonal
    d = diag(A) # n vector
    full_fact_space += n
    full_fact_fevals += n
    error = norm(d,1)
    piv = [1:n] # n vector
    full_fact_space += n
    diag_blocks = zeros(num_blocks)
    j = 1
    while error > tol
        for i=j:num_blocks
            upper_i = min(piv[i]*block_size,n)
            diag_blocks[i] = sum(d[(piv[i]-1)*block_size+1:upper_i])
        end
        idx = indmax(diag_blocks) + j - 1
        piv[j],piv[idx] = piv[idx],piv[j]
        lower_j = (piv[j]-1)*block_size+1
        upper_j = min(piv[j]*block_size,n)
        range_j = lower_j:upper_j
        block_size_j = upper_j - lower_j + 1
        push!(L,zeros(n,block_size_j)) # n by block_size_j matrix
        full_fact_space += n*block_size_j
        for i=j:num_blocks
            lower_i = (piv[i]-1)*block_size+1
            upper_i = min(piv[i]*block_size,n)
            range_i = lower_i:upper_i
            block_size_i = upper_i - lower_i + 1
            # evaluate block
            S = A[range_i,range_j]
            full_fact_fevals += block_size_i*block_size_j
            for k=1:j-1
                S -=
                   L[k][range_i,1:block_size]*
                   L[k][range_j,1:block_size]'
            end
            if i==j
                G = chol(S,:L)
            else
                G = S/L[j][range_j,1:block_size_j]'
                d[range_i] -= [(G[m,:]*G[m,:]')[1] for m=1:block_size_i]
            end
            L[j][range_i,1:block_size_j] = G
        end
        error = 0
        for i=j+1:num_blocks
            upper_i = min(piv[i]*block_size,n)
            error += sum(d[(piv[i]-1)*block_size+1:upper_i])
        end
        j += 1
    end
    L2 = zeros(n,min(block_size*(j-1),n))
    full_fact_space += n*block_size*(j-1)
    for p=1:j-1
        upper_p = min(p*block_size,n)
        lower_p = (p-1)*block_size+1
        block_size_p = upper_p - lower_p + 1
        L2[:,lower_p:upper_p] = L[p]
    end
    return L2, piv, full_fact_space, full_fact_fevals
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

function sorted_blocked_full_fact(A,tol,block_size)
    # Computes the sorted blocked full factorization of A
    full_fact_space = 0
    full_fact_fevals = 0
    n = size(A,1)
    num_blocks = convert(Int64,div(n,block_size)) + (convert(Bool,mod(n,block_size)) ? 1 : 0)
    L = Array{Float64,2}[]
    # evaluate diagonal
    d = diag(A) # n vector
    full_fact_space += n
    full_fact_fevals += n
    error = norm(d,1)
    piv = [1:n] # n vector
    full_fact_space += n
    diag_blocks = zeros(num_blocks)
    j = 1
    while error > tol
        for i=j:num_blocks
            upper_i = min(piv[i]*block_size,n)
            diag_blocks[i] = sum(d[(piv[i]-1)*block_size+1:upper_i])
        end
        idx = indmax(diag_blocks) + j - 1
        piv[j],piv[idx] = piv[idx],piv[j]
        lower_j = (piv[j]-1)*block_size+1
        upper_j = min(piv[j]*block_size,n)
        range_j = lower_j:upper_j
        block_size_j = upper_j - lower_j + 1
        push!(L,zeros(n,block_size_j)) # n by block_size_j matrix
        full_fact_space += n*block_size_j
        for i=j:num_blocks
            lower_i = (piv[i]-1)*block_size+1
            upper_i = min(piv[i]*block_size,n)
            range_i = lower_i:upper_i
            block_size_i = upper_i - lower_i + 1
            # evaluate block
            S = A[range_i,range_j]
            full_fact_fevals += block_size_i*block_size_j
            for k=1:j-1
                S -=
                   L[k][range_i,1:block_size]*
                   L[k][range_j,1:block_size]'
            end
            if i==j
                G = chol(S,:L)
            else
                G = S/L[j][range_j,1:block_size_j]'
                d[range_i] -= [(G[m,:]*G[m,:]')[1] for m=1:block_size_i]
            end
            L[j][range_i,1:block_size_j] = G
        end
        error = 0
        for i=j+1:num_blocks
            upper_i = min(piv[i]*block_size,n)
            error += sum(d[(piv[i]-1)*block_size+1:upper_i])
        end
        j += 1
    end
    L2 = zeros(n,min(block_size*(j-1),n))
    full_fact_space += n*block_size*(j-1)
    for p=1:j-1
        upper_p = min(p*block_size,n)
        lower_p = (p-1)*block_size+1
        block_size_p = upper_p - lower_p + 1
        L2[:,lower_p:upper_p] = L[p]
    end
    return L2, piv, full_fact_space, full_fact_fevals
end


function blocked_struct_fact(A,tol,block_size,nb)
    # Computes the blocked full factorization of A
    fact_space = 0
    fact_fevals = 0
    n = size(A,1)
    nb = convert(Int64,sqrt(n))
    nsym = convert(Int64,nb*(nb+1)/2) # size(A,1)
    nsym_div_block_size = convert(Int64,div(nsym,block_size))
    extra_block = convert(Bool,mod(nsym,block_size)) ? 1 : 0
    num_blocks = nsym_div_block_size+extra_block
    L = Array{Float64,2}[]
    # evaluate diagonal
    d = struct_diag(A) # nsym vector
    fact_space += nsym
    fact_fevals += nsym
    error = norm(d,1)
    piv = [1:nsym] # nsym vector
    fact_space += nsym
    diag_blocks = zeros(num_blocks)
    fact_space += num_blocks
    j = 1
    while error > tol
        for i=j:num_blocks
            upper_i = min(piv[i]*block_size,nsym)
#            println("upper_i: ",upper_i)
#            println("nsym: ",nsym)
#            println("length(d): ",length(d))
            diag_blocks[i] = sum(d[(piv[i]-1)*block_size+1:upper_i])
        end
        idx = indmax(diag_blocks) + j - 1
        piv[j],piv[idx] = piv[idx],piv[j]
        lower_j = (piv[j]-1)*block_size+1
        upper_j = min(piv[j]*block_size,nsym)
        range_j = lower_j:upper_j
        block_size_j = upper_j - lower_j + 1
        push!(L,zeros(nsym,block_size_j)) # n by block_size_j matrix
        fact_space += nsym*block_size_j
        for i=j:num_blocks
            println("i: ",i,", j: ",j)
            lower_i = (piv[i]-1)*block_size+1
            upper_i = min(piv[i]*block_size,nsym)
            range_i = lower_i:upper_i
            block_size_i = upper_i - lower_i + 1
            # evaluate block
            #S = A[range_i,range_j]
            S = struct_block(A,lower_i,lower_j,block_size_i,block_size_j)
            fact_fevals += block_size_i*block_size_j
            for k=1:j-1
                block_size_k = min(piv[k]*block_size,nsym)-(piv[k]-1)*block_size
                S -=
                   L[k][range_i,1:block_size_k]*
                   L[k][range_j,1:block_size_k]'
            end
            if i==j
                println("rank(S): ",rank(S))
                G = chol(S,:L)
            else
                G = S/L[j][range_j,1:block_size_j]'
                d[range_i] -= [(G[m,:]*G[m,:]')[1] for m=1:block_size_i]
            end
            L[j][range_i,1:block_size_j] = G
        end
        error = 0
        for i=j+1:num_blocks
            upper_i = min(piv[i]*block_size,nsym)
            error += sum(d[(piv[i]-1)*block_size+1:upper_i])
        end
        println("Error: ",error)
        j += 1
    end
    L2 = zeros(nsym,min(block_size*(j-1),nsym))
    fact_space += nsym*block_size*(j-1)
    for p=1:j-1
        upper_p = min(p*block_size,nsym)
        lower_p = (p-1)*block_size+1
        block_size_p = upper_p - lower_p + 1
        L2[:,lower_p:upper_p] = L[p]
    end
    return L2, piv, fact_space, fact_fevals
end

function sorted_blocked_struct_fact(A,tol,block_size,nb)
    # Computes the blocked full factorization of A
    fact_space = 0
    fact_fevals = 0
    n = size(A,1)
    nb = convert(Int64,sqrt(n))
    nsym = convert(Int64,nb*(nb+1)/2) # size(A,1)
    nsym_div_block_size = convert(Int64,div(nsym,block_size))
    extra_block = convert(Bool,mod(nsym,block_size)) ? 1 : 0
    num_blocks = nsym_div_block_size+extra_block
    L = Array{Float64,2}[]
    # evaluate diagonal
    d = struct_diag(A) # nsym vector
    fact_space += nsym
    fact_fevals += nsym
    error = norm(d,1)
    piv = [1:nsym] # nsym vector
    fact_space += nsym
    diag_blocks = zeros(num_blocks)
    fact_space += num_blocks
    j = 1
    while error > tol
        for i=j:num_blocks
            upper_i = min(piv[i]*block_size,nsym)
#            println("upper_i: ",upper_i)
#            println("nsym: ",nsym)
#            println("length(d): ",length(d))
            diag_blocks[i] = sum(d[(piv[i]-1)*block_size+1:upper_i])
        end
        idx = indmax(diag_blocks) + j - 1
        piv[j],piv[idx] = piv[idx],piv[j]
        lower_j = (piv[j]-1)*block_size+1
        upper_j = min(piv[j]*block_size,nsym)
        range_j = lower_j:upper_j
        block_size_j = upper_j - lower_j + 1
        push!(L,zeros(nsym,block_size_j)) # n by block_size_j matrix
        fact_space += nsym*block_size_j
        for i=j:num_blocks
            println("i: ",i,", j: ",j)
            lower_i = (piv[i]-1)*block_size+1
            upper_i = min(piv[i]*block_size,nsym)
            range_i = lower_i:upper_i
            block_size_i = upper_i - lower_i + 1
            # evaluate block
            #S = A[range_i,range_j]
            S = struct_block(A,lower_i,lower_j,block_size_i,block_size_j)
            fact_fevals += block_size_i*block_size_j
            for k=1:j-1
                block_size_k = min(piv[k]*block_size,nsym)-(piv[k]-1)*block_size
                S -=
                   L[k][range_i,1:block_size_k]*
                   L[k][range_j,1:block_size_k]'
            end
            if i==j
                println("rank(S): ",rank(S))
                G = chol(S,:L)
            else
                G = S/L[j][range_j,1:block_size_j]'
                d[range_i] -= [(G[m,:]*G[m,:]')[1] for m=1:block_size_i]
            end
            L[j][range_i,1:block_size_j] = G
        end
        error = 0
        for i=j+1:num_blocks
            upper_i = min(piv[i]*block_size,nsym)
            error += sum(d[(piv[i]-1)*block_size+1:upper_i])
        end
        println("Error: ",error)
        j += 1
    end
    L2 = zeros(nsym,min(block_size*(j-1),nsym))
    fact_space += nsym*block_size*(j-1)
    for p=1:j-1
        upper_p = min(p*block_size,nsym)
        lower_p = (p-1)*block_size+1
        block_size_p = upper_p - lower_p + 1
        L2[:,lower_p:upper_p] = L[p]
    end
    return L2, piv, fact_space, fact_fevals
end


function test_facts()
    tol = 1e-6
    nbs = [34, 48, 68]
    block_sizes = [96, 256, 768]
    ranks = [200, 300, 400]
    for i=1:length(nbs)#nb=nbs
        nb = nbs[i]
        block_size = block_sizes[i]
        rank = ranks[i]
        println("nb: ",nb)
        #    n = convert(Int64,nb*(nb+1)/2)
        #    A = randn(n,n)
        #    A = A*A'
        A = RandERIMatrix(nb,rank)

#        println("Blocked full fact")
#        tic()
#        L2, piv, full_fact_space, full_fact_fevals = blocked_full_fact(A,tol,block_size)
#        toc()
#        error = norm(A-L2*L2',1)
#        println("Error: ", error)

        println("Blocked struct fact")
        tic()
        L2, piv, full_fact_space, full_fact_fevals = blocked_struct_fact(A,tol,block_size,nb)
        toc()
        J = QsymTransform(L2,1:convert(Int64,nb*(nb+1)/2),nb)
        error = norm(A-J*J',1)
        println("Error: ", error)

        println("Full fact")
        tic()
        L2, piv, full_fact_space, full_fact_fevals = full_fact(A,tol)
        toc()
        error = norm(A[piv,piv]-L2*L2',1)
        println("Error: ", error)

        println("Struct fact")
        tic()
        L, piv, struct_fact_space, struct_fact_fevals = struct_fact(A,tol)
        toc()
        J = QsymTransform(L,piv,nb)
        error = norm(A - J*J',1)
        println("Error: ", error)

        println("LAPACK full fact")
        Ain = deepcopy(A)
        tic()
        L, piv, struct_fact_space, struct_fact_fevals = lapack_full_fact(Ain,tol)
        toc()
        error = norm(A[piv,piv]-L*L',1)
        println("Error: ", error)

        println("LAPACK struct fact")    
        n2,n2 = size(A)
        n = convert(Int64,sqrt(n2))
        Q = Qnn(n)
        nsym = convert(Int64,n*(n+1)/2)
        Qsym = Q[:,1:nsym]
        Asym = Qsym'*A*Qsym
        tic()
        L2, piv2, struct_fact_space, struct_fact_fevals = lapack_struct_fact(Asym,tol)
        toc()
        J = QsymTransform(L2,piv2,nb)
        error = norm(A - J*J',1)
        println("Error: ", error)    

        C = randn(nb,nb)

        println("Full transformation")
        tic()
        B, full_trans_space = full_trans(L,piv,C)
        full_trans_time = toc()
        full_trans_error = norm(kron(C,C)*A*kron(C,C)' - B*B',1)
        println("error: ",full_trans_error)

        println("Structured transformation")
        tic()
        B, struct_trans_space = struct_trans(L2,piv2,C)
        struct_trans_time = toc()
        struct_trans_error = norm(kron(C,C)*A*kron(C,C)' - B*B',1)
        println("error: ",struct_trans_error)
    end    
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

function testBlockDiag()
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
end

function TestStructTransform()
    n = 3
    r = 3
    A = RandERIMatrix(n,r)
    Q = Qnn(n)
    nsym = convert(Int64,n*(n+1)/2)
    Qsym = Q[:,1:nsym]
    Asym = Qsym'*A*Qsym
    C = randn(n,n)
    Csym = Q'*kron(C,C)*Q
    println("Csym: ")
    display(Csym)
end

function dpstrfL(A,tol,nb)
    # Compute the Cholesky factorization P' * A * P = L * L'
    n = size(A,1)
    piv = [1:n]
    work = zeros(2n)
    dstop = tol

    i__2 = n
    i__1 = nb
    for k=1:nb:n
                                                              
        # Account for last block not being NB wide */

        # Computing MIN */
        i__3 = nb
        i__4 = n - k + 1
        jb = min(i__3,i__4)

        # Set relevant part of first half of WORK to zero, */
        # holds dot products */

        i__3 = n
        for i__=k:i__3
	    work[i__] = 0
            # L150: */
        end
	i__3 = k + jb - 1
        for j=k:i__3            

            # Find pivot, test for exit, else swap rows and columns */
            # Update dot products, compute possible pivots which are */
            # stored in the second half of WORK */
            
            i__4 = n
            for i__=j:i__4
                
	        if (j > k) 
                    # Computing 2nd power */
	            d__1 = a[i__,j-1]
	            work[i__] += d__1 * d__1
	        end		                                      
	        work[n + i__] = a[i__,i__] - work[i__]
                
                # L160: */
            end
	    if (j > 1) 
	        maxlocvar = (n << 1) - (n + j) + 1
	        itemp = LAPACK.dmaxloc_(work[n + j], maxlocvar)
	        pvt = itemp + j - 1
	        ajj = work[n + pvt]
	        if (ajj <= dstop || LAPACK.disnan_(&ajj)) 
		    a[j + j * a_dim1] = ajj;
                    # break
		    #goto L190;
	        end	
	    end	    
            
	    if (j != pvt) 
                
                # Pivot OK, so can now swap pivot rows and columns */
                
	        a[pvt,pvt] = a[j,j]
	        i__4 = j - 1;
	        dswap_(&i__4, &a[j + a_dim1], lda, &a[pvt + a_dim1], 
		       lda);
	        if (pvt < n) 
		    i__4 = n - pvt;
		    dswap_(&i__4, &a[pvt + 1 + j * a_dim1], &c__1, 
                           &a[pvt + 1 + pvt * a_dim1], &c__1);
                end
	        i__4 = pvt - j - 1;
	        dswap_(&i__4, &a[j + 1 + j * a_dim1], &c__1, 
                       &a[pvt+(j + 1) * a_dim1], lda);
                
                # Swap dot products and PIV */

	        dtemp = work[j];
	        work[j] = work[pvt];
	        work[pvt] = dtemp;
	        itemp = piv[pvt];
	        piv[pvt] = piv[j];
	        piv[j] = itemp;
            end
	    ajj = sqrt(ajj);
	    a[j + j * a_dim1] = ajj;

            # Compute elements J+1:N of column J. */

	    if (j < n) 
	        i__4 = n - j;
	        i__5 = j - k;
	        dgemv_("No Trans", &i__4, &i__5, &c_b22, 
                       &a[j + 1 + k * a_dim1], lda, &a[j + k * a_dim1],
                       lda, &c_b24, &a[j + 1 + j * a_dim1], &c__1);
	        i__4 = n - j;
	        d__1 = 1. / ajj;
	        dscal_(&i__4, &d__1, &a[j + 1 + j * a_dim1], &c__1);
	    end
            # L170: */
        end

        # Update trailing matrix, J already incremented */

        if (k + jb <= n)
	    i__3 = n - j + 1;
	    dsyrk_("Lower", "No Trans", &i__3, &jb, &c_b22, 
                   &a[j + k * a_dim1], lda, 
                   &c_b24, &a[j + j * a_dim1], lda);
        end
        # L180: */
    end
    return L,piv
end

function lapack_full_fact(A,tol)
#    Ain = deepcopy(A)
    Aout, piv, rank, info = LAPACK.pstrf!('L', A, tol)
    L = tril(Aout)
    L = L[:,1:rank]
    return L, piv, 0, 0
end

function lapack_struct_fact(Asym,tol)
    Aout, piv, rank, info = LAPACK.pstrf!('L', Asym, tol)
    L = tril(Aout)
    L = L[:,1:rank]
    p
    return L, piv, 0, 0
end

function TestDsptrfL()
    n = 9
    nb = 3
    A = randn(n,6)
    A = A*A'
    L,piv = dpstrfL(A,tol,nb)
    error =  norm(A[piv,piv]-L*L')
    println("Error: ",error)
end

#TestIntegralTransformation()
#TestStructTransform()
#TestDsptrfL()
test_facts()