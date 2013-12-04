function S = StructReshapeQskewTimesVector(v,n)
% Returns Qskew*v computed without matrix multiplication reshaped as n x n
% S is skew-symmetric since Qskew is an orthogonal basis for reshaped
% skewmetric matricies

S = zeros(n);
% The lower triangular part of S is enumerated as 
% the elements of v
S(find(tril(ones(n),-1))) = v/sqrt(2);
S = S - S';