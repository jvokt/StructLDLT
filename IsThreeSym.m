function b = IsThreeSym(A)
% Returns true if A is three symmetric (has symmetric blocks, 
% block symmetry, and perfect shuffle permutation symmetry)

[n2,~] = size(A);
n = sqrt(n2);
% Testing symmetric blocks
for k=1:n
    for l=1:n
        I = n*(k-1)+1:n*k;
        J = n*(l-1)+1:n*l;
        block = A(I,J);
        sym_block_err = norm(block - block');
    end
end
% Testing block symmetry
for k=1:n
    for l=k+1:n
        I = n*(k-1)+1:n*k;
        J = n*(l-1)+1:n*l;
        block1 = A(I,J);
        block2 = A(J,I);
        block_symm_err = norm(block1 - block2);
    end
end
% Testing perfect shuffle permutation symmetry
perf_shuff_err = norm(A-A(PerfShuff(n,n),PerfShuff(n,n)));

b = sym_block_err + block_symm_err + perf_shuff_err < 1e-10;