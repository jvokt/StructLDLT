function QTA = Qsym(Asym,n,~)
% Returns Qsym*Asym without matrix multiplication
% Software dependencies: StructReshapeQsymTimesVector
[~,n2] = size(Asym);
QTA = zeros(n^2,n2);
% Computes the resulting matrix column by column
for k=1:n2
    QTA(:,k) = reshape(StructReshapeQsymTimesVector(Asym(:,k),n),n^2,1);
end