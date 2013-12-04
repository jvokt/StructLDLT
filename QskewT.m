function QTA = QskewT(A,n,nskew)
% Returns Qskew'*A without matrix multiplication
% Software dependencies: QskewTransposedTimesVector
[~,n2] = size(A);
QTA = zeros(nskew,n2);
% Computes the resulting matrix column by column
for k=1:n2
    QTA(:,k) = QskewTransposedTimesVector(A(:,k),n,nskew);
end