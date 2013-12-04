function w = QskewTransposedTimesVector(v,n,~)
% Returns Qsym'*v without matrix multiplication

Ak = reshape(v,n,n);
% The first n elements of w are the diagonal elements of Ak
%w(1:n) = diag(Ak);
AkT = Ak';
% v1 is a vector of the the lower triangular elements of Ak
v1 = Ak(find(tril(ones(n),-1)));
% v2 is a vector of the the upper triangular elements of Ak
v2 = AkT(find(tril(ones(n),-1)));
% w is sqrt(2) times the average of v1 and v2
w = (v1 - v2)/sqrt(2);