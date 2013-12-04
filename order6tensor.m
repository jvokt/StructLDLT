n = 3;
k = 1;
T1 = zeros(n,n,n,n,n,n);
for j3=1:n
   for j2=1:n
      for j1=1:n
          for i3=1:n
             for i2=1:n
                for i1=1:n
                    T1(i1,i2,i3,j1,j2,j3) = k;
                    k = k+1;
                end
            end
          end
      end
   end
end
k = 1;
T2 = zeros(n,n,n,n,n,n);
for i3=1:n
   for i2=1:n
      for i1=1:n
         for j3=1:n
            for j2=1:n
                for j1=1:n
                    T2(i1,i2,i3,j1,j2,j3) = k;
                    k = k+1;
                end
            end
          end
      end
   end
end

A1 = TenToMat(T1,1:3,4:6);
A1sub = A1(1:18,1:18);
A2 = TenToMat(T1,4:6,1:3);
%A2(1:18,1:18)
error = norm(A1'-A2)
%A1v = Vec(A1);
%A2v = Vec(A2);
%p = PerfShuff(3^3,3^3);
%norm(A1v-A2v(p));

k = 1;
T3 = zeros(n,n,n,n,n,n);
for j3=1:n
   for j2=1:n
      for j1=1:n
         for i3=1:n
            for i1=1:n
               for i2=1:n
                   T3(i1,i2,i3,j1,j2,j3) = k;
                   k = k+1;
               end
            end
         end
      end
   end
end

A3 = TenToMat(T3,1:3,4:6);
A3sub = A3(1:18,1:18);
I3 = eye(3,3);
P3 = eye(9,9);
P3 = P3(:,PerfShuff(3,3));
IP3 = kron(I3,P3);
error = norm(IP3*A1 - A3)

k = 1;
T4 = zeros(n,n,n,n,n,n);
for j3=1:n
   for j2=1:n
      for j1=1:n
         for i2=1:n
            for i3=1:n
               for i1=1:n
                   T4(i1,i2,i3,j1,j2,j3) = k;
                   k = k+1;
               end
            end
         end
      end
   end
end

A4 = TenToMat(T4,1:3,4:6);
A4sub = A4(1:18,1:18)
PI3 = kron(P3,I3);
error = norm(PI3*A1 - A4)
