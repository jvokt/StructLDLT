function QTA = QsymT(A,n,nsym)
% Returns Qsym'*A without matrix multiplication
% Software dependencies: QsymTransposedTimesVector
[~,n2] = size(A);
QTA = zeros(nsym,n2);
% Computes the resulting matrix column by column
for k=1:n2
    QTA(:,k) = QsymTransposedTimesVector(A(:,k),n,nsym);
end