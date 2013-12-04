function S = StructReshapeQsymTimesVector(v,n)
% Returns Qsym*v computed without matrix multiplication reshaped as n x n
% S is symmetric since Qsym is an orthogonal basis for reshaped
% symmetric matricies

% The diagonal elements of S are the first n elements of v
D = diag(v(1:n));
S = zeros(n);
% The lower triangular part of S is enumerated as 
% the remaining elements of v
S(find(tril(ones(n),-1))) = v(n+1:end)/sqrt(2);
S = S + D + S';