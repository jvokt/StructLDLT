function QTA = Qskew(Askew,n,~)
% Returns Qskew*A without matrix multiplication
% Software dependencies: StructReshapeQskewTimesVector
[~,n2] = size(Askew);
QTA = zeros(n^2,n2);
% Computes the resulting matrix column by column
for k=1:n2
    QTA(:,k) = reshape(StructReshapeQskewTimesVector(Askew(:,k),n),n^2,1);
end