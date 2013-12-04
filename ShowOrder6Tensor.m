function ShowOrder6Tensor()

n = 3;
n3 = n^3;
A = rand(n3);
A = A*A';
A = Order6symmetrize(A,n);

[Qe,Qa,Qs] = Qorder6(n);

Q = [Qe, Qa, Qs];

%Q(abs(Q)<1e-10) = 0;
%spy(Q)
%pause
B = Q'*A*Q;
%As = Qs'*A*Qs;
%[U,T] = schur(As);
%[Y I] = sort(diag(T));
%T = T(I,I);
%U(abs(U)<1e-6) = 0;
%U = U(:,I)
%spy(U)

%norm(As(1:8,1:8)-As(9:16,9:16))
%norm(As(1:8,1:8)-As(1:8,1:8)')
%norm(As(1:8,9:16)+As(1:8,9:16)')
%norm(As(1:8,9:16)-As(9:16,1:8)')
%size(As)
%B = reshape(As,4,4,4,4);
%B = B + permute(B, [2; 1; 4; 3]);
%B = reshape(B,16,16)

rs = length(As);

%norm(As(rs/2+1:end,1:rs/2)+As(rs/2+1:end,1:rs/2)')

Qg = [eye(rs/2,rs/2), eye(rs/2,rs/2); eye(rs/2,rs/2), -eye(rs/2,rs/ ...
                                                  2)];


Qg'*As*Qg



%{
ii = zeros(rs,1);
jj = zeros(rs,1);
ss = zeros(rs,1);
nz = 1;
c = 1;
for i=1:rs/2
    ii(nz) = i;
    jj(nz) = c;
    ss(nz) = 1/sqrt(2);
    nz = nz+1;
    
    ii(nz) = i+rs/2;
    jj(nz) = c;
    ss(nz) = 1/sqrt(2);
    nz = nz+1;
    
    c = c+1;
end

Qf = sparse(ii,jj,ss,rs,rs/2)

ii = zeros(rs,1);
jj = zeros(rs,1);
ss = zeros(rs,1);
nz = 1;
c = 1;
for i=1:rs/2
    ii(nz) = i;
    jj(nz) = c;
    ss(nz) = 1/sqrt(2);
    nz = nz+1;
    
    ii(nz) = i+rs/2;
    jj(nz) = c;
    ss(nz) = -1/sqrt(2);
    nz = nz+1;
    
    c = c+1;
end

Qg = sparse(ii,jj,ss,rs,rs/2)

Qd = [Qf, Qg];

Ad = Qd'*As*Qd;

Ad(abs(Ad)<1e-10) = 0;

spy(Ad)
%}
%[V,D] = eig(As);
%evalnorm = norm(T-D)
%evecnorm = norm(U-V)
%B = Q'*A*Q;
%B(abs(B)<1e-10) = 0;
%spy(B)
%diag(As)
%Qs(:,[1 9])
%rank(U)
%diag(T)
%U(:,1:6)

%As

count = 0;

for i=1:rs
    for j=1:i
        for k=1:rs
            for l=1:k
                if (i~=k && j~=l ) && abs(As(i,j) - As(k,l)) < 1e-10 && ...
                       As(i,j) ~= 0
                    %fprintf('(%d,%d)==(%d,%d)\n',i,j,k,l);
                    count = count + 1;
                end
            end
        end
    end
end

%disp(count)

%diag(A)

%{
for i = [2 4 10]
    for j = [2 4 10]
        fprintf('(%d,%d) is %d\n',i,j,A(i,j));
    end
end


As11 = A(2,2)*(2/3) + A(2,4)*sqrt(2/3)*(-1/sqrt(6)) + A(2,10)* ...
       sqrt(2/3)*(-1/sqrt(6)) + A(4,2)*(-1/sqrt(6))*sqrt(2/3) + ...
       A(4,4)*(1/6) + A(4,10)*(1/6) + A(10,2)*(-1/sqrt(6))*sqrt(2/3) ...
       + A(10,4)*(1/6) + A(10,10)*(1/6);

err1 = abs(As11-As(1,1))

As99 = A(4,4)*(1/2) + A(4,10)*(-1/2) + A(10,4)*(-1/2) + A(10,10)*(1/2);

err9 = abs(As99-As(9,9))

err19 = abs(As11-As99)
%}
function As = Order6symmetrize(A,n)
% Permutations on three indices
G = [1, 1, 2, 2, 3, 3;
     2, 3, 1, 3, 1, 2;
     3, 2, 3, 1, 2, 1];

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

function [Qe,Qa,Qs] = Qorder6(n)
% Q is the n^3 x n^3 orthonormal matrix which block diagonalizes the 
% tensor unfolding
% Q has (1/3) * n * (11 * n^2 - 9*n + 1) non-zero entries
% Qe is dimension n + n*(n-1) + n*(n-1)*(n-2)/6
% Qa is dimension n*(n-1)*(n-2)/6
% Qs is dimension n^3 - n - n*(n-1) - n*(n-1)*(n-2)/3

c = 1;
nz = 1;
nnze = n+n*(n-1)*3+n*(n-1)*(n-2);
ii = zeros(nnze,1);
jj = zeros(nnze,1);
ss = zeros(nnze,1);

% Specifies n non-zeros
for i=1:n
    ii(nz) = col([i i i],[n n n]);
    jj(nz) = c;
    ss(nz) = 1;
    nz = nz + 1;
    c = c + 1;
end

% Specifies n*(n-1)*3 non-zeros
for j=1:n
    for i=[1:j-1,j+1:n]
        ii(nz) = col([i j j],[n n n]);
        jj(nz) = c;
        ss(nz) = 1/sqrt(3);
        nz = nz + 1;
        ii(nz) = col([j i j],[n n n]);
        jj(nz) = c;
        ss(nz) = 1/sqrt(3);
        nz = nz + 1;
        ii(nz) = col([j j i],[n n n]);
        jj(nz) = c;
        ss(nz) = 1/sqrt(3);
        nz = nz + 1;
        c = c + 1;
    end
end

% Specifies n*(n-1)*(n-2) non-zeros
for k=1:n
    for j=1:k-1
        for i=1:j-1
            ii(nz) = col([i j k],[n n n]);
            jj(nz) = c;
            ss(nz) = 1/sqrt(6);
            nz = nz + 1;
            ii(nz) = col([i k j],[n n n]);
            jj(nz) = c;
            ss(nz) = 1/sqrt(6);
            nz = nz + 1;
            ii(nz) = col([j i k],[n n n]);
            jj(nz) = c;
            ss(nz) = 1/sqrt(6);
            nz = nz + 1;
            ii(nz) = col([j k i],[n n n]);
            jj(nz) = c;
            ss(nz) = 1/sqrt(6);
            nz = nz + 1;
            ii(nz) = col([k i j],[n n n]);
            jj(nz) = c;
            ss(nz) = 1/sqrt(6);
            nz = nz + 1;
            ii(nz) = col([k j i],[n n n]);
            jj(nz) = c;
            ss(nz) = 1/sqrt(6);
            nz = nz + 1;
            c = c + 1;
        end
    end
end

re = n + n*(n-1) + n*(n-1)*(n-2)/6;%nchoosek(n,1) + 2*nchoosek(n,2) + nchoosek(n,3);

Qe = sparse(ii,jj,ss,n^3,re);

nnza = n*(n-1)*(n-2);
ii = zeros(nnza,1);
jj = zeros(nnza,1);
ss = zeros(nnza,1);

nz = 1;
c = 1;

% Specifies n*(n-1)*(n-2) non-zeros
for k=1:n
    for j=1:k-1
        for i=1:j-1
            ii(nz) = col([i j k],[n n n]);
            jj(nz) = c;
            ss(nz) = 1/sqrt(6);
            nz = nz + 1;
            ii(nz) = col([i k j],[n n n]);
            jj(nz) = c;
            ss(nz) = -1/sqrt(6);
            nz = nz + 1;
            ii(nz) = col([j i k],[n n n]);
            jj(nz) = c;
            ss(nz) = -1/sqrt(6);
            nz = nz + 1;
            ii(nz) = col([j k i],[n n n]);
            jj(nz) = c;
            ss(nz) = 1/sqrt(6);
            nz = nz + 1;
            ii(nz) = col([k i j],[n n n]);
            jj(nz) = c;
            ss(nz) = 1/sqrt(6);
            nz = nz + 1;
            ii(nz) = col([k j i],[n n n]);
            jj(nz) = c;
            ss(nz) = -1/sqrt(6);
            nz = nz + 1;
            c = c + 1;
        end
    end
end

ra = n*(n-1)*(n-2)/6;%nchoosek(n,3);
Qa = sparse(ii,jj,ss,n^3,ra);

nnzs = (5/3)*(n-1)*n*(n+1);
ii = zeros(nnzs,1);
jj = zeros(nnzs,1);
ss = zeros(nnzs,1);

rs = n^3 - re - ra;

nz = 1;
c = 1;

% Specifies (5/3)*(n-1)*n*(n+1) non-zeros
for k=1:n-1
    for j=k:n
        for i=k+1:n
            ii(nz) = col([i j k],[n n n]);
            jj(nz) = c;
            ss(nz) = sqrt(2/3);
            nz = nz + 1;
            ii(nz) = col([j k i],[n n n]);
            jj(nz) = c;
            ss(nz) = -1/sqrt(6);
            nz = nz + 1;
            ii(nz) = col([k i j],[n n n]);
            jj(nz) = c;
            ss(nz) = -1/sqrt(6);
            nz = nz + 1;

            ii(nz) = col([j k i],[n n n]);
            jj(nz) = c+rs/2;
            ss(nz) = 1/sqrt(2);
            nz = nz + 1;
            ii(nz) = col([k i j],[n n n]);
            jj(nz) = c+rs/2;
            ss(nz) = -1/sqrt(2);
            nz = nz + 1;
            c = c + 1;
        end
    end
end

Qs = sparse(ii,jj,ss,n^3,rs);
end
end