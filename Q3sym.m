function Q = Q3sym(n)
% Q is the n^2xn^2 orthonormal matrix defined in the UCLA pdf

Q = sparse(n^2,n^2);
for i=1:n
    Q(1+(n+1)*(i-1),i) = 1;
end

k = n+1;
p = n*(n+1)/2+1;
for i=1:n
    for j=i+1:n
        if (i==j)
            Q(1+(n+1)*(i-1),i) = 1;
        else
            Q(i+(j-1)*n,k) = 1/sqrt(2);
            Q(j+(i-1)*n,k) = 1/sqrt(2);
            Q(i+(j-1)*n,p) = -1/sqrt(2);
            Q(j+(i-1)*n,p) = 1/sqrt(2);
            k = k+1;
            p = p+1;
        end
    end
end

