function [Qe,Qa,Qs] = Qorder6(n)
% Q is the n^3 x n^3 orthonormal matrix defined which block diagonalizes the 
% tensor unfolding
% Q has (1/3) * n * (11 * n^2 - 9*n + 1) non-zero entries

c = 1;
nz = 1;

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

re = nchoosek(n,1) + 2*nchoosek(n,2) + nchoosek(n,3);

Qe = sparse(ii,jj,ss,n^3,re);

ii = [];
jj = [];
ss = [];

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

ra = nchoosek(n,3);
Qa = sparse(ii,jj,ss,n^3,ra);

ii = [];
jj = [];
ss = [];

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