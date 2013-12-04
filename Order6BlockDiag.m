function [A, Ae, Aa, As,Pe,Pa,Ps,Qe,Qa,Qs] = Order6BlockDiag(n)
% Block diagonalization of a reshaped order-6 ERI tensor

% Permutations on three indices
G = [1, 1, 2, 2, 3, 3;
     2, 3, 1, 3, 1, 2;
     3, 2, 3, 1, 2, 1];

% Characters for the three representations
chi = [1,  1,  1,  1,  1,  1;
       1, -1, -1,  1,  1, -1;
       2,  0,  0, -1, -1,  0];

% Form a ERI3-symmetric matrix
A = rand(n^3);
A = (A+A')/2;
A = symmetrize3(A,n);

sym_err = zeros(6,1);

for i=1:6
    idx = [G(:,i),G(:,i)+3];
    T = reshape(A,[n n n n n n]);
    sym_err(i) = norm(A-reshape(permute(T,idx),n^3,n^3));
end

% Form projectors
Pe = form_P(1,n);
Pa = form_P(2,n);
Ps = form_P(3,n);

spy(Ps)
Ps(:,2)
Ps(:,4)
Ps(:,10)

B = Ps([2 4 10],[2 4 10])
[Q R] = rrqr(B,'tol',1e-15);

% % Check commutation relations
% fprintf('-----\n')
% fprintf('Comm Pe: %e\n', norm(Pe*As-As*Pe)/norm(As));
% fprintf('Comm Pa: %e\n', norm(Pa*As-As*Pa)/norm(As));
% fprintf('Comm Ps: %e\n', norm(Ps*As-As*Ps)/norm(As));
% 
% % Check that they're orthogonal projectors
% se = svd(Pe);
% sa = svd(Pa);
% ss = svd(Ps);
% fprintf('-----\n')
% fprintf('Orth proj Pe?  %e\n', sum(abs(se.*(1-se))));
% fprintf('Orth proj Pa?  %e\n', sum(abs(sa.*(1-sa))));
% fprintf('Orth proj Ps?  %e\n', sum(abs(ss.*(1-ss))));

% Look at the ranks
fprintf('-----\n')
fprintf('Rank Pe: %g\n', rank(Pe))
fprintf('Rank Pa: %g\n', rank(Pa))
fprintf('Rank Ps: %g\n', rank(Ps))

% Check that they have different ranges
% fprintf('-----\n')
% fprintf('Mutual orth: %e %e %e\n', norm(Pe*Pa), norm(Pe*Ps), norm(Pa*Ps));

% Super-symmetrize A
function As = symmetrize3(A,n)
  A = reshape(A, n,n,n, n,n,n);
  As = A;
  idx1 = [1;2;3];
  idx2 = [4;5;6];
  for k = 2:6
    idx = [idx1(G(:,k)); idx2(G(:,k))];
    As = As + permute(A, idx);
  end
  As = reshape(As, n^3, n^3)/6;
end

% Apply the representation rep to a order-3 tensor V in dimension n
function PV = apply_P(rep,n,V)
  V = reshape(V, n,n,n, numel(V)/n^3);
  PV = chi(rep,1)*V;
  for k = 2:6
    PV = PV + chi(rep,k)*permute(V, [G(:,k); 4]);
  end
  %PV = PV*(chi(rep,1)/6);
  PV = PV/6;
  PV = reshape(PV, n^3, numel(V)/n^3);
end

% Form the projector associated with a representation
function P = form_P(rep,n)
  P = apply_P(rep, n, eye(n^3));
  P = reshape(P, n^3, n^3);
end

% 
% d = 6;
% 
% n3 = n^3;
% A = rand(n3);
% A = A+A';
% 
% % Permutations on three indices
% %G = perms(1:d/2)
% G = [1, 1, 2, 2, 3, 3;
%      2, 3, 1, 3, 1, 2;
%      3, 2, 3, 1, 2, 1];
% 
% B = reshape(1:n3,n,n,n);
% 
% P = zeros(n3,factorial(d/2));
% 
% for g=1:length(G)
%     P(:,g) = Vec(permute(B,G(g,:)));
% end
% 
% % Symmetrize over group
% B = zeros(n3);
% for g = 1:length(G)
%   B = B + A( P(:,g), P(:,g) );
% end
% A = B;
% 
% size(A)
% 
% %A = B/6;
% %B(:,1)'
% %T = MatToTen(B, [1 2 3], [4 5 6],[3 3 3 3 3 3])
% 
% % Characters for the three representations
% chi = [1,  1,  1,  1,  1,  1;
%        1, -1, -1,  1,  1, -1;
%        2,  0,  0, -1, -1,  0];
% 
% % Form the three canonical projectors
% Pe = zeros(n3);
% Pa = zeros(n3);
% Ps = zeros(n3);
% I  = eye(n3);
% for g = 1:length(G)
%   Pe = Pe + chi(1,g)*I(P(:,g),:); 
%   Pa = Pa + chi(2,g)*I(P(:,g),:); 
%   Ps = Ps + chi(3,g)*I(P(:,g),:); 
% end
% Pe = Pe/6;
% Pa = Pa/6;
% Ps = Ps/3;
% 
% % Sanity check
% 
% fprintf('\nSingular values of canonical projectors\n');
% se = svd(Pe);
% sa = svd(Pa);
% ss = svd(Ps);
% re = rank(Pe,1e-10)
% ra = rank(Pa,1e-10)
% rs = rank(Ps,1e-10)
% 
% reas = rank([Pe Pa Ps])
% rea = rank([Pe Pa])
% res = rank([Pe Ps])
% ras = rank([Pa Ps])
% 
% %fprintf('\nDimension of range space of projectors\n');
% sum(se);
% sum(sa);
% sum(ss);
% 
% %fprintf('\nCheck invariance\n');
% norm(Pa'*A*Pe)/norm(A);
% norm(Ps'*A*Pe)/norm(A);
% norm(Ps'*A*Pa)/norm(A);
% 
% fprintf('\nClassify eigenvectors\n');
% [V,D] = eig(A);
% [sum((Pe*V).^2); sum((Pa*V).^2); sum((Ps*V).^2); diag(D)'];

%[Qe,Re] = rrqr(Pe,'rnk',5);
[Qe,Re] = rrqr(Pe,'tol',1e-10);
%[~,r] = size(Qe);
%spy(Qe)

%[Qa,Ra] = rrqr(Pa,'rnk',4);
[Qa,Ra] = rrqr(Pa,'tol',1e-10);
%[~,r] = size(Qa);
%spy(Qa)

%[Qs,Rs] = rrqr(Ps,'rnk',18);
[Qs,Rs] = rrqr(Ps,'tol',1e-10);
%[~,r] = size(Qs);

%fprintf('Mutual orth: %e %e %e\n', norm(Qe'*Qa), norm(Qe'*Qs), norm(Qa'*Qs));

Ae = Qe'*A*Qe;
Aa = Qa'*A*Qa;
As = Qs'*A*Qs;
end