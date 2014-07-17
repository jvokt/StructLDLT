using PyPlot

function tuned_algorithm_comparison()
    tol = 1e-6
    molecules = ["HF","NH3","H2O2","N2H4","C2H5OH"]
    num_basis_fns = [34, 48, 68, 82, 123]
    ranks = [200, 300, 400, 500, 750]
    block_sizes = [64, 96, 96, 96, 64]
    nmols = length(molecules)
    algs = ["Full tensor transformation",
            "Tuned Algorithm 4.1 plus Cholesky factor transform"]
    style = ["wo","go"]
    full_timings = zeros(nmols)
    chol_timings = zeros(nmols)
    qchol_timings = zeros(nmols)
    full_errors = zeros(nmols)
    tuned_errors = zeros(nmols)
    for i=1:nmols
        nbf = num_basis_fns[i]
        n = nbf^2
        r = ranks[i]
        nb = block_sizes[i]
        println("n: ",n)
        A = randn(n,r)
        A = A*A'
        C = randn(nbf,nbf)
        tic()
        M = unfact_trans(A,C)
        full_timings[i] = toc()
#        full_errors[i] = trace(kron(C,C)*A*kron(C,C)'-M)
#        println("Error: ",full_errors[i])
	         
        tic()
        L = blocked_full_fact_row_major(A,tol,nb)
        G = half_trans(L,C)
        chol_timings[i] = toc()
	Asym = randn(int(n*(n+1)/2),r)
	Asym = Asym*Asym'
        tic()
        Lsym = blocked_full_fact_row_major(Asym,tol,nb)
        G = half_trans(L,C)
        qchol_timings[i] = toc()
#        tuned_errors[i] = trace(kron(C,C)*A*kron(C,C)'-G*G')
#        println("Error: ", tuned_errors[i])
    end
    if false   
    figure()
    plot(num_basis_fns[1:nmols],full_timings,style[1],label=algs[1])
    plot(num_basis_fns[1:nmols],tuned_timings,style[2],label=algs[2])
    legend()
    xlabel("n (number of basis functions)")
    ylabel("Time in seconds")
    title("State-of-the-art versus new tuned algorithm timing")
    savefig(string("State-of-the-art versus new tuned algorithm timing",".png"))
    figure()
    semilogy(num_basis_fns[1:nmols],full_errors,style[1],label=algs[1])
    semilogy(num_basis_fns[1:nmols],tuned_errors,style[2],label=algs[2])
    legend()
    xlabel("n (number of basis functions)")
    ylabel("Trace norm error")
    title("State-of-the-art versus new tuned algorithm error")
    savefig(string("State-of-the-art versus new tuned algorithm error",".png"))
    end
    println(full_timings)
    println(tuned_timings)
end

function brmllrrchol(A,tol,block_size)
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

function tune_memory_layout()
    tol = 1e-6
    molecules = ["HF","NH3","H2O2","N2H4","C2H5OH"]
    num_basis_fns = [34, 48, 68, 82, 123]
    ranks = [200, 300, 400, 500, 750]
    block_sizes = [32, 64, 96, 128, 192, 256, 384, 512]
    nmols = length(molecules)
    facts = ["Algorithm 4.1, diagonal blocks",
             "Algorithm 4.1",
             "Algorithm 4.1, diagonal blocks, row-major block order",
             "Algorithm 4.1, row-major block order"]
    style = ["ro","bo","go","wo"]
    num_facts = length(facts)
    timings = zeros(nmols,length(block_sizes),num_facts)
    for i=1:nmols
        nbf = num_basis_fns[i]
        n = nbf^2
        r = ranks[i]
        println("n: ",n)
        A = randn(n,r)
        A = A*A'
        for b = 1:length(block_sizes)
            nb = block_sizes[b]
            println("block size: ",nb)

            # Blocked pivoted Cholesky (precompute diagonal blocks)
            tic()
            L2, piv, full_fact_space, full_fact_fevals = blocked_full_fact(A,tol,nb)
            timings[i,b,1] = toc()

            # Blocked pivoted Cholesky (precompute diagonal)
            tic()
            L2, piv, full_fact_space, full_fact_fevals = blocked_full_fact_diag(A,tol,nb)
            timings[i,b,2] = toc()

            # Blocked pivoted Cholesky (precompute diagonal blocks, row-major)
            tic()
            L2, piv, full_fact_space, full_fact_fevals = blocked_full_fact_row_major(A,tol,nb)
            timings[i,b,3] = toc()

            # Blocked pivoted Cholesky (precompute diagonal, row-major)
            tic()
            L2, piv, full_fact_space, full_fact_fevals = blocked_full_fact_diag_row_major(A,tol,nb)
            timings[i,b,4] = toc()
        end
    end
    
    for i=1:nmols
        figure()
        for f=1:num_facts
            plot(block_sizes,vec(timings[i,:,f]),style[f],label=facts[f])
        end
        legend()
        xlabel("Block size")
        ylabel("Time in seconds")
        title(string("All timings for ", molecules[i]))
        savefig(string("All timings for ", molecules[i],".png"))
    end
end

function TestIntegralTransformation()

    molecules = ["HF","NH3","H2O2","N2H4","C2H5OH"]
    num_basis_fns = [34, 48, 68, 82, 123]
    ranks = [200, 300, 400, 500, 750]
    nmols = 2#length(num_basis_fns)
    gen_plots = false
    tol = 1e-6
    
    unfact_trans_time = zeros(nmols)
    unfact_trans_space = zeros(nmols)
    unfact_trans_error = zeros(nmols)
    unfact_trans_fevals = zeros(nmols)

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
        unfact_trans_time[i] = toc()
        unfact_error[i] = trace(kron(C,C)*A*kron(C,C)'-M)
        println("error: ",unfact_error[i])

        # Full factorization
        println("Full factorization")
        tic()
        L, piv, full_fact_space[i], full_fact_fevals[i] = 
        full_fact(A,tol)
        full_fact_time[i] = toc()
        full_fact_error[i] = trace(A[piv,piv] - L*L')
        println("error: ", full_fact_error[i])

        # Full transformation
        println("Full transformation")
        tic()
        B, full_trans_space[i] = full_trans(L,piv,C)
        full_trans_time[i] = toc()
        full_trans_error[i] = trace(kron(C,C)*A*kron(C,C)' - B*B')
        println("error: ",full_trans_error[i])

        println("Total time for full algorithm: ",
                full_fact_time[i]+full_trans_time[i])
        
        # Structured factorization
        println("Structured factorization")
        tic()
        L, piv, struct_fact_space[i], struct_fact_fevals[i] = 
        struct_fact(A,tol)
        struct_fact_time[i] = toc()
        J = QsymTransform(L,piv,n)
        struct_fact_error[i] = trace(A - J*J')
        println("error: ",struct_fact_error[i])
        
        # Structured transformation
        println("Structured transformation")
        tic()
        B, struct_trans_space[i] = struct_trans(L,piv,C)
        struct_trans_time[i] = toc()
        struct_trans_error[i] = trace(kron(C,C)*A*kron(C,C)' - B*B')
        println("error: ",struct_trans_error[i])

        println("Total time for structured algorithm: ",
                struct_fact_time[i]+struct_trans_time[i])
    end
    
    if gen_plots
        plotComparison(num_basis_fns[1:nmols],
                       full_fact_time,
                       "Full factorization",
                       struct_fact_time,
                       "Structured factorization",
                       "Number of basis functions",
                       "Time (seconds)",
                       "Time to factorize")
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
                       "Error (1-norm)",
                       "Factorization Error")
        plotComparison(num_basis_fns[1:nmols],
                       full_fact_fevals,
                       "Full factorization",
                       struct_fact_fevals,
                       "Structured factorization",
                       "Number of basis functions",
                       "Fevals",
                       "Fevals to factorize")
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
                       "Space (in Float64s)",
                       "Space to transform")
        plotComparison(num_basis_fns[1:nmols],
                       full_trans_error,
                       "Full transformation",
                       struct_trans_error,
                       "Structured transformation",
                       "Number of basis functions",
                       "Error (trace norm)",
                       "Transformation Error")
        plotComparison3(num_basis_fns[1:nmols],
                       full_fact_time+full_trans_time,
                       "Full factorization plus transformation",
                       struct_fact_time+struct_trans_time,
                       "Structured factorization plus transformation",
                       unfact_trans_time,
                       "Unfactorized transformation",
                       "Number of basis functions",
                       "Time (seconds)",
                       "Time to factorize and transform")
        plotComparison3(num_basis_fns[1:nmols],
                       full_fact_space+full_trans_space,
                       "Full factorization plus transformation",
                       struct_fact_space+struct_trans_space,
                       "Structured factorization plus transformation",
                       unfact_trans_space,
                       "Unfactorized transformation",
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

function plotComparison3(x,y1,l1,y2,l2,y3,l3,xlab,ylab,titl)
    figure()
    plot(x,y1,"ro",label=l1)
    plot(x,y2,"bo",label=l2)
    plot(x,y3,"go",label=l3)
    legend()
    xlabel(xlab)
    ylabel(ylab)
    title(titl)
    savefig(string(titl,".png"))
end

function plotComparison4(x,y1,l1,y2,l2,y3,l3,y4,l4,xlab,ylab,titl)
    figure()
    plot(x,y1,"ro",label=l1)
    plot(x,y2,"bo",label=l2)
    plot(x,y3,"go",label=l3)
    plot(x,y4,"wo",label=l4)
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

function unblocked_psi4(A,tol)
    # Computes the full factorization of A
    n = size(A,1)
    col = zeros(n) # n vector
    L = Array{Float64,1}[]
    m = 1
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
    L2 = zeros(n,m-1)
    for p=1:m-1
        L2[:,p] = L[p]
    end
    return L2
end

function blocked_full_fact(A,tol,block_size)
    n = size(A,1)
    num_blocks = convert(Int64,div(n,block_size)) + (convert(Bool,mod(n,block_size)) ? 1 : 0)
    D = Array{Float64,2}[]
    for i=1:num_blocks
        range_i = (i-1)*block_size+1:min(i*block_size,n)
        push!(D,A[range_i,range_i])
    end
    L = Array{Float64,2}[]
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
        push!(L,zeros(n,block_size_j))
        G = lapack_full_fact(D[piv[j]],tol)
        r = size(G,2)
        L[j][range_j,1:r] = G
        for i=j+1:num_blocks
            lower_i = (piv[i]-1)*block_size+1
            upper_i = min(piv[i]*block_size,n)
            range_i = lower_i:upper_i
            block_size_i = upper_i-lower_i+1
            # GEMM update of subdiagonal blocks
            S = A[range_i,range_j]
            for k=1:j-1
                S -=
                L[k][range_i,1:block_size]*
                L[k][range_j,1:block_size]'
            end
            G = S/L[j][range_j,1:r]'
            # SYRK update of diagonal blocks
            D[piv[i]] -= G*G'
            L[j][range_i,1:r] = G
        end
        error = sum(map(trace,D[piv[j+1:num_blocks]]))
        j+=1
    end
    L2 = zeros(n,min(block_size*(j-1),n))
    for p=1:j-1
        upper_p = min(p*block_size,n)
        lower_p = (p-1)*block_size+1
        L2[:,lower_p:upper_p] = L[p]
    end
    return L2,piv,0,0
end

function blocked_full_fact_row_major(A,tol,block_size)
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

function blocked_full_fact_diag(A,tol,block_size)
    n = size(A,1)
    num_blocks = convert(Int64,div(n,block_size)) + (convert(Bool,mod(n,block_size)) ? 1 : 0)
    d = reshape(vcat(diag(A),zeros(num_blocks*block_size-n)),block_size,num_blocks)    
    L = Array{Float64,2}[]
    piv = [1:num_blocks]
    error = sum(d)
    j = 1
    while error > tol
        idx = indmax(sum(d[:,piv[j:num_blocks]],1))+j-1
        piv[j],piv[idx] = piv[idx],piv[j]
        lower_j = (piv[j]-1)*block_size+1
        upper_j = min(piv[j]*block_size,n)
        range_j = lower_j:upper_j
        block_size_j = upper_j-lower_j+1
        push!(L,zeros(n,block_size_j))
        r = block_size
        for i=j:num_blocks
            lower_i = (piv[i]-1)*block_size+1
            upper_i = min(piv[i]*block_size,n)
            range_i = lower_i:upper_i
            block_size_i = upper_i-lower_i+1
            S = A[range_i,range_j]
            for k=1:j-1
                S -=
                L[k][range_i,1:block_size]*
                L[k][range_j,1:block_size]'
            end
            if i==j
                G = lapack_full_fact(S,tol)
                r = size(G,2)
            else
                G = S/L[j][range_j,1:r]'
                d[1:block_size_i,piv[i]] -= 
                [(G[m,:]*G[m,:]')[1] for m=1:block_size_i]
            end
            L[j][range_i,1:r] = G
        end
        error = sum(d[:,piv[j:num_blocks]])
        j+=1
    end
    L2 = zeros(n,min(block_size*(j-1),n))
    for p=1:j-1
        lower_p = (p-1)*block_size+1
        L2[:,lower_p:lower_p+size(L[p],2)-1] = L[p]
    end
    return L2,piv,0,0
end

function blocked_full_fact_diag_row_major(A,tol,block_size)
    n = size(A,1)
    num_blocks = convert(Int64,div(n,block_size)) + (convert(Bool,mod(n,block_size)) ? 1 : 0)
    d = reshape(vcat(diag(A),zeros(num_blocks*block_size-n)),block_size,num_blocks)    
    L = [Array{Float64,2}[] for i=1:num_blocks]
    piv = [1:num_blocks]
    error = sum(d)
    j = 1
    while error > tol
        idx = indmax(sum(d[:,piv[j:num_blocks]],1))+j-1
        piv[j],piv[idx] = piv[idx],piv[j]
        lower_j = (piv[j]-1)*block_size+1
        upper_j = min(piv[j]*block_size,n)
        range_j = lower_j:upper_j
        block_size_j = upper_j-lower_j+1
        r = block_size
        for i=j:num_blocks
            lower_i = (piv[i]-1)*block_size+1
            upper_i = min(piv[i]*block_size,n)
            range_i = lower_i:upper_i
            block_size_i = upper_i-lower_i+1
            S = A[range_i,range_j]
            for k=1:j-1
                S -= L[piv[i]][k]*L[piv[j]][k]'
            end
            if i==j
                G = lapack_full_fact(S,tol)
                r = size(G,2)
            else
                G = S/L[piv[j]][j]'
                d[1:block_size_i,piv[i]] -= 
                [(G[m,:]*G[m,:]')[1] for m=1:block_size_i]
            end
            push!(L[piv[i]],G)
        end
        error = sum(d[:,piv[j:num_blocks]])
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
    return L2,piv,0,0
end

function blocked_struct_fact(A,tol,block_size,nb)
    n = convert(Int64,nb*(nb+1)/2)
    num_blocks = convert(Int64,div(n,block_size)) + (convert(Bool,mod(n,block_size)) ? 1 : 0)
    D = Array{Float64,2}[]
    for i=1:num_blocks
        lower_i=(i-1)*block_size+1
        block_size_i=min(i*block_size,n)-lower_i+1
        push!(D,struct_block(A,lower_i,lower_i,block_size_i,block_size_i))
#        range_i = (i-1)*block_size+1:min(i*block_size,n)
#        push!(D,A[range_i,range_i])#zeros(block_size,block_size))
#        D[i] = A[range_i,range_i]
    end
    L = Array{Float64,2}[]
    piv = [1:num_blocks]
    piv2 = [1:block_size]
    error = sum(map(trace,D))
    j = 1
    while error > tol
        idx = indmax(map(trace,D[piv[j:num_blocks]]))+j-1
        piv[j],piv[idx] = piv[idx],piv[j]
        lower_j = (piv[j]-1)*block_size+1
        upper_j = min(piv[j]*block_size,n)
        range_j = lower_j:upper_j
        block_size_j = upper_j-lower_j+1
        push!(L,zeros(n,block_size_j))
        G = lapack_full_fact(D[piv[j]],1e-12)
#        println(size(G))
#        println(size(range_j))
        r = size(G,2)
        L[j][range_j,1:r] = G
        for i=j+1:num_blocks
            lower_i = (piv[i]-1)*block_size+1
            upper_i = min(piv[i]*block_size,n)
            range_i = lower_i:upper_i
            block_size_i = upper_i-lower_i+1
#            S = A[range_i,range_j]
            S = struct_block(A,lower_i,lower_j,block_size_i,block_size_j)
            for k=1:j-1
                S -=
                L[k][range_i,1:block_size]*
                L[k][range_j,1:block_size]'
            end
            G = S/L[j][range_j,1:r]'
            D[piv[i]] -= G*G'
            L[j][range_i,1:r] = G
        end
        error = sum(map(trace,D[piv[j+1:num_blocks]]))
        j+=1
    end
    L2 = zeros(n,min(block_size*(j-1),n))
    for p=1:j-1
        upper_p = min(p*block_size,n)
        lower_p = (p-1)*block_size+1
        block_size_p = upper_p - lower_p + 1
        L2[:,lower_p:upper_p] = L[p]
    end
    return L2,piv,0,0
end

function blocked_full_fact3(A,tol,block_size)
    # Computes the blocked full factorization of A
    full_fact_space = 0
    full_fact_fevals = 0
    n = size(A,1)
    num_blocks = convert(Int64,div(n,block_size)) + (convert(Bool,mod(n,block_size)) ? 1 : 0)
    L = Array{Float64,2}[]
    # evaluate diagonal
#    d = diag(A) # n vector
    D = Array{Float64,2}[]
    for i=1:num_blocks
        upper_i = min((i-1)*block_size,n)
        block_size_i = upper_i - (i-1)*block_size
        push!(D,zeros(block_size_i,block_size_i)) # n by block_size_j matrix
        D[i] = A[(i-1)*block_size+1:upper_i,(i-1)*block_size+1:upper_i]
    end
    full_fact_space += n
    full_fact_fevals += n
    error = sum(map(trace,D))
    piv = [1:num_blocks] # num_blocks vector
    full_fact_space += n
    block_trace = zeros(num_blocks)
    j = 1
    while error > tol
#        for i=j:num_blocks
#            block_trace[i] = trace(D[piv[i]])
#        end
#        idx = indmax(block_trace[j:num_blocks]) + j - 1
        idx = indmax(map(trace,D[piv[j:num_blocks]])) + j - 1
        piv[j],piv[idx] = piv[idx],piv[j]
        lower_j = (piv[j]-1)*block_size+1
        upper_j = min(piv[j]*block_size,n)
        range_j = lower_j:upper_j
        block_size_j = upper_j - lower_j + 1
        push!(L,zeros(n,block_size_j)) # n by block_size_j matrix
        full_fact_space += n*block_size_j
        for i=j:num_blocks
            # evaluate block
#            println("Rank of S: ",rank(S))
            if last_block_column && i==j+1
                lower_i = (piv[i]-1)*block_size+1-r
                upper_i = min(piv[i]*block_size,n)
                range_i = lower_i:upper_i
                block_size_i = upper_i - lower_i + 1
                lower_j = (piv[j]-1)*block_size+1
                upper_j = min(piv[j]*block_size,n)-r
                range_j = lower_j:upper_j
                block_size_j = upper_j - lower_j + 1
            else
                lower_i = (piv[i]-1)*block_size+1
                upper_i = min(piv[i]*block_size,n)
                range_i = lower_i:upper_i
                block_size_i = upper_i - lower_i + 1
            end

            if i==j
                r = rank(D[piv[i]])
                G = chol(D[piv[i]][1:r,1:r],:L)
                if r < block_size
                    last_block_column = true
                end
            else
                S = A[range_i,range_j]
                full_fact_fevals += block_size_i*block_size_j
                for k=1:j-1
                    S -=
                    L[k][range_i,1:block_size]*
                    L[k][range_j,1:block_size]'
                end
                G = S/L[j][range_j,1:block_size_j]'
                D[piv[i]] -= G*G'
            end
            L[j][range_i,1:block_size_j] = G
        end
        error = sum(map(trace,D[piv[j+1:num_blocks]]))
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

function symmat(A)                                       
    n2 = size(A,1)
    n = convert(Int64,sqrt(n2))
    nsym = convert(Int64,n*(n+1)/2)
    Q = Qnn(n)
    Qsym = Q[:,1:nsym]
    if A == A'
        return struct_block(A,1,1,nsym,nsym)
    else
        return Qsym'*A*Qsym
    end
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

function blocked_full_fact2(A,tol,block_size)
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
#        p = sortperm(d,rev=true)
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

function blocked_struct_fact2(A,tol,block_size,nb)
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
#            println("i: ",i,", j: ",j)
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

function half_trans(L,C)
    # Computes the factorized full transformation of A
    n,n = size(C)
    n2,r = size(L)
    B1 = C*reshape(L,n,n*r)
    B2 = C*reshape(permutedims(reshape(B1,n,n,r),[2 1 3]),n,n*r)
    B3 = reshape(permutedims(reshape(B2,n,n,r),[2 1 3]),n2,r)
    return B3
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
    # d = diag(Asym)
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
        # Asym[piv[1:nsym],piv[m]]
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
        G[:,p] = L[p]#[piv]
    end
    return G#(G,piv)
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
    Aout, piv, rank, info = LAPACK.pstrf!('L', A, tol)
    L = tril(Aout)
    L[piv,1:rank] = L[:,1:rank]
    return L[:,1:rank]
end

function lapack_struct_fact(Asym,tol)
    Aout, piv, rank, info = LAPACK.pstrf!('L', Asym, tol)
    L = tril(Aout)
    L = L[:,1:rank]
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

function randInRange(n,num)
    indices = zeros(num); x = 0;
    for i=1:num
        x = floor(n*rand())+1;
        while true
            if !(x in indices)
                break
            else
                x = round(n*rand())
            end
        end
        indices[i] = x
    end
    return indices
end

function test_blocked_pivoted_cholesky()
    n = 34^2
    r = 34^2
    tol = 1e-6
    A = randn(n,r)
    A = A*A'
    block_size = 32
    num_blocks = convert(Int64,div(n,block_size)) + (convert(Bool,mod(n,block_size)) ? 1 : 0)
    println("Block size: ",block_size)
    println("Number of blocks: ",num_blocks)
    D = Array{Float64,2}[]
    for i=1:num_blocks
        push!(D,zeros(block_size,block_size))
        range_i = (i-1)*block_size+1:min(i*block_size,n)
        D[i] = A[range_i,range_i]
#        println("size(D[i])",size(D[i]))
    end
    L = Array{Float64,2}[]
    piv = [1:num_blocks]
    piv2 = [1:block_size]
    error = sum(map(trace,D))
    j = 1
    while error > tol
        idx = indmax(map(trace,D[piv[j:num_blocks]]))+j-1
        piv[j],piv[idx] = piv[idx],piv[j]
        lower_j = (piv[j]-1)*block_size+1
        upper_j = min(piv[j]*block_size,n)
        range_j = lower_j:upper_j
        block_size_j = upper_j-lower_j+1
        push!(L,zeros(n,block_size_j))
#        println("size(D[piv[j]])",size(D[piv[j]]))
        G = LazyCholesky(D[piv[j]],1e-16)
        r = size(G,2)
        L[j][range_j,1:r] = G
        for i=j+1:num_blocks
            lower_i = (piv[i]-1)*block_size+1
            upper_i = min(piv[i]*block_size,n)
            range_i = lower_i:upper_i
            block_size_i = upper_i-lower_i+1
            S = A[range_i,range_j]
            for k=1:j-1
                S -=
                L[k][range_i,1:block_size]*
                L[k][range_j,1:block_size]'
            end
            G = S/L[j][range_j,1:r]'
            D[piv[i]] -= G*G'
            L[j][range_i,1:r] = G
        end
        error = sum(map(trace,D[piv[j+1:num_blocks]]))
        j+=1
    end
    L2 = zeros(n,min(block_size*(j-1),n))
    for p=1:j-1
        upper_p = min(p*block_size,n)
        lower_p = (p-1)*block_size+1
        block_size_p = upper_p - lower_p + 1
        L2[:,lower_p:upper_p] = L[p]
    end
    
    error = trace(A-L2*L2')
    println("\nError: ",error)
#    println("num_blocks in piv[1:j-1]?: ",num_blocks in piv[1:j-1])
end

function test_facts()
    tol = 1e-6
    nbs = [34, 48, 68]
    block_sizes = [96, 256, 768]
    ranks = [200, 300, 400]
    for i=1:length(nbs)
        nb = nbs[i]
        block_size = block_sizes[i]
        rank = ranks[i]
        println("nb: ",nb)
        A = RandERIMatrix(nb,rank)

        println("Blocked full fact")
        tic()
        L2, piv, full_fact_space, full_fact_fevals = blocked_full_fact(A,tol,block_size)
        toc()
        error = trace(A-L2*L2')
        println("Error: ", error)

        println("Blocked full fact (block diag)")
        tic()
        L2, piv, full_fact_space, full_fact_fevals = blocked_full_fact_block_diag(A,tol,block_size)
        toc()
        error = trace(A-L2*L2')
        println("Error: ", error)

        println("Blocked struct fact")
        tic()
        L2, piv, full_fact_space, full_fact_fevals = blocked_struct_fact(A,tol,block_size,nb)
        toc()
        J = QsymTransform(L2,1:convert(Int64,nb*(nb+1)/2),nb)
        error = norm(A-J*J',1)
#        error = norm(A-L2*L2',1)
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

function tune_blocked_facts()
    tol = 1e-6
    molecules = ["HF","NH3","H2O2","N2H4","C2H5OH"]
    num_basis_fns = [34, 48, 68, 82, 123]
    ranks = [200, 300, 400, 500, 750]
    block_sizes = [32, 64, 96, 128, 192, 256, 384, 512]
    nmols = 2#length(molecules)
#    facts = ["Blocked pivoted Cholesky (diagonal blocks, column-major)",
    facts = ["Algorithm 4.1 with precomputing diagonal blocks",
#             "Blocked pivoted Cholesky (diagonal, column-major)",
             "Algorithm 4.1",
             "LAPACK pstrf",
             "Unfactorized transformation",
             "Blocked pivoted Cholesky (diagonal blocks, row-major)",
             "Blocked pivoted Cholesky (diagonal, row-major)"]
    style = ["ro","bo","go","wo","r+","b+"]
    num_facts = length(facts)
    timings = zeros(nmols,length(block_sizes),num_facts)
    for i=1:nmols
        nbf = num_basis_fns[i]
        n = nbf^2
        r = ranks[i]
        C = randn(nbf,nbf)
        println("n: ",n)
        A = randn(n,r)
        A = A*A'
        for b = 1:length(block_sizes)
            nb = block_sizes[b]
            println("block size: ",nb)

            # Blocked pivoted Cholesky (precompute diagonal blocks)
            println(facts[1])
            tic()
            L2, piv, full_fact_space, full_fact_fevals = blocked_full_fact(A,tol,nb)
            timings[i,b,1] = toc()
            error = trace(A-L2*L2')
            println("Error: ", error)

            # Blocked pivoted Cholesky (precompute diagonal)
            println(facts[2])
            tic()
            L2, piv, full_fact_space, full_fact_fevals = blocked_full_fact_diag(A,tol,nb)
            timings[i,b,2] = toc()
            error = trace(A-L2*L2')
            println("Error: ", error)

            # LAPACK pstrf
            println(facts[3])
 #           Ain = deepcopy(A)
            tic()
 #           L2 = lapack_full_fact(Ain,tol)
            timings[i,b,3] = toc()
 #           error = trace(A-L2*L2')
 #           println("Error: ", error)                         

            # Unfactorized transformation
            println(facts[4])
            tic()
#            M = unfact_trans(A,C)
            timings[i,b,4] = toc()
#            error = trace(kron(C,C)*A*kron(C,C)'-M)
#            println("Error: ", error)

            # Blocked pivoted Cholesky (precompute diagonal blocks, row-major)
            println(facts[5])
            tic()
 #           L2, piv, full_fact_space, full_fact_fevals = blocked_full_fact_row_major(A,tol,nb)
            timings[i,b,5] = toc()
#            println("size(A[piv,piv]): ",size(A[piv,piv]))
#            println("size(A[piv,piv]): ",size(L2*L2'))
#            error = trace(A-L2*L2')
#            println("Error: ", error)

            # Blocked pivoted Cholesky (precompute diagonal, row-major)
            println(facts[6])
            tic()
#            L2, piv, full_fact_space, full_fact_fevals = blocked_full_fact_diag_row_major(A,tol,nb)
            timings[i,b,6] = toc()
#            error = trace(A-L2*L2')
#            println("Error: ", error)
        end
    end
    
    for i=1:nmols
        figure()
#        for f=1:num_facts
        for f=[1 2]
            if f == 3 || f == 4
                plot(block_sizes,ones(length(block_sizes))*timings[i,1,f],style[f],label=facts[f])
            else
                plot(block_sizes,vec(timings[i,:,f]),style[f],label=facts[f])
            end                
        end
        legend()
        xlabel("Block size")
        ylabel("Time in seconds")
        title(string("Algorithm 4.1 with precomputing diagonal blocks timings for ", molecules[i]))
        savefig(string("Algorithm 4.1 with precomputing diagonal blockstimings for ", molecules[i],".png"))
    end
end

function test_blocked_full_fact_row_major()
    nbf = 34
    nb = 32
    n = nbf^2
    r = 200
    tol = 1e-6
    C = randn(nbf,nbf)
    println("n: ",n)
    A = randn(n,r)
    A = A*A'
    tic()
    L2, piv, full_fact_space, full_fact_fevals = blocked_full_fact_row_major(A,tol,nb)
    toc()
    error = trace(A-L2*L2')
    println("Error: ", error)
end

function profile_facts()
    tol = 1e-6
    molecules = ["HF","NH3","H2O2","N2H4","C2H5OH"]
    num_basis_fns = [34, 48, 68, 82, 123]
    ranks = [200, 300, 400, 500, 750]
    nmols = 2
    routines = ["Unblocked Psi4 decomposition",
                "LAPACK PSTRF decomposition",
                "Cholesky factor transformation",
                "Full tensor transformation"]
    style = ["ro","bo","go","wo"]
    num_routines = length(routines)
    timings = zeros(nmols,num_routines)
    for i=1:nmols
        nbf = num_basis_fns[i]
        n = nbf^2
        r = ranks[i]
        C = randn(nbf,nbf)
        println("n: ",n)
        A = randn(n,r)
        A = A*A'

        # Unblocked Psi4 decomposition
        println(routines[1])
        tic()
        L = unblocked_psi4(A,tol)
        timings[i,1] = toc()
        error = trace(A-L*L')
        println("Error: ", error)

        # LAPACK PSTRF decomposition
        println(routines[2])
        Ain = deepcopy(A)
        tic()
        L = lapack_full_fact(Ain,tol)
        timings[i,2] = toc()
        error = trace(A-L*L')
        println("Error: ", error)                         

        # Cholesky factor transformation
        println(routines[3])
        tic()
        G = half_trans(L,C)
        timings[i,3] = toc()
        error = trace(kron(C,C)*A*kron(C,C)'-G*G')
        println("Error: ", error)

        # Full tensor transformation
        println(routines[4])
        tic()
        B = unfact_trans(A,C)
        timings[i,4] = toc()
        error = trace(kron(C,C)*A*kron(C,C)'-B)
        println("Error: ", error)
    end
    
    figure()
    for r=1:num_routines
        plot(num_basis_fns[1:nmols],timings[:,r],style[r],label=routines[r])
    end
    legend()
    xlabel("Number of basis functions (n)")
    ylabel("Time in seconds")
    title("Profiling unblocked Psi4 decomposition")
    savefig("Profiling unblocked Psi4 decomposition.png")
end

function gen_mol_data(molecule)
    # Generate "molecule".dat from Psi4
    
end

function parse_data_gen_ten(molecule)
    # Parse "molecule".dat for two-electron integrals, generate "molecule".ten
end

function read_ten_unfold_plot(molecule)
    # Read "molecule".ten into a tensor, unfold as a matrix, plot eigenvalues
end

TestIntegralTransformation()
#TestStructTransform()
#TestDsptrfL()
#test_facts()
#test_blocked_pivoted_cholesky()
#tune_blocked_facts()
#test_blocked_full_fact_row_major()
#profile_facts()
#tune_memory_layout()
#tuned_algorithm()