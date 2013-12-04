function [B, sig, C, tau] = StructLDLTOrder6(A,tol)
% Factorizes A according to the Order 6 tensor analogy
% Software dependencies: LazyLDLT, QsymT, StructReshapeQsymTimesVector,
% QskewT, StructReshapeQskewTimesVector

[n3,~] = size(A);
n = n3^(1/3);

[Qe,Qa,Qs] = Qorder6(n);

% Computing Asym = Qsym'*A*Qsym
% QsymTA = Qsym'*A
% QsymTA = QsymT(A,n,nsym);
% Asym = (Qsym'*(Qsym'*A)')'
% Asym = QsymT(QsymTA',n,nsym)';
% Lazy-evaluation LDLT factorization
[Lsym,sig,psym] = LazyLDLT(@F, tol, nsym, Asym);
% rsym is the revealed rank of Asym
rsym = length(sig);
v = zeros(rsym,1);
B = cell(rsym);
for k=1:rsym
    % v = P'*Lsym(:,k)
    v(psym) = Lsym(:,k);
    % Reshapes Qsym*v computed without matrix multiplication
    B{k} = StructReshapeQsymTimesVector(v,n);
end
% Computing Askew = Qskew'*A*Qskew
% QskewTA = Qskew'*A
QskewTA = QskewT(A,n,nskew);
% Askew = (Qskew'*(Qskew'*A)')'
Askew = QskewT(QskewTA',n,nskew)';
% Lazy-evaluation LDLT factorization
[Lskew,tau,pskew] = LazyLDLT(@F, tol, nskew, Askew);
% rskew is the revealed rank of Askew
rskew = length(tau);
v = zeros(rskew,1);
C = cell(rskew);
for k=1:rskew
    % v = P'*Lskew(:,k)
    v(pskew) = Lskew(:,k);
    % Reshapes Qskew*v computed without matrix multiplication
   C{k} = StructReshapeQskewTimesVector(v,n);
end
