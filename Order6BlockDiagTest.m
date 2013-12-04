n = 3;
[A, Ae, Aa, As,Pe,Pa,Ps,Qe,Qa,Qs] = Order6BlockDiag(n);
Qs
% Qe(abs(Qe)<1e-10) = 0;
% Qa(abs(Qa)<1e-10) = 0;
% Qs(abs(Qs)<1e-10) = 0;
% spy(Qe);
% rank(Qe)
% figure
% spy(Qa);
% rank(Qa)
% figure
% spy(Qs);
% rank(Qs)

%[Qe2,Qa2,Qs2] = Qorder6(n);
%Qe2(abs(Qe2)<1e-10) = 0;
%Qa2(abs(Qa2)<1e-10) = 0;
%Qs2(abs(Qs2)<1e-10) = 0;
% figure
% spy(Qe2);
% rank(full(Qe2))
% figure
% spy(Qa2);
% rank(full(Qa2))
% figure
% spy(Qs2);
% rank(full(Qs2))

% Q = [Qe, Qa, Qs];
%Q2 = [Qe2, Qa2, Qs2];
% I = Q'*Q;
% I(abs(I)<1e-10) = 0;
% I2 = Q2'*Q2;
% I2(abs(I2)<1e-10) = 0;
% spy(I)
% figure
% spy(I2)
% figure
% spy(Qs)
% figure
% spy(Qs2)
% B = Q'*A*Q;
%size(Q2)
%size(A)
%B2 = Q2'*A*Q2;
% B(B<1e-10) = 0;
%B2(B2<1e-10) = 0;
% figure
% spy(B)
%figure
%spy(B2)
% for j=1:rank(Qs)
%     j
%     x = Qs(:,j);
%     x(x<1e-15) = 0;
%     x = find(x);
%     inv_col(x(1),[n n n])
% end