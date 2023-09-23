function W = chorin_conv_matrix(c4n,n4e,Grads_T,Vol_T,u)
[nC,d] = size(c4n); nE = size(n4e,1);
ctr = 0; ctr_max = (d+1)^2*nE;
I = zeros(ctr_max,1); J = zeros(ctr_max,1); 
m_loc = (ones(d+1,d+1)+eye(d+1))/((d+1)*(d+2));
X_W = zeros(ctr_max,1); 
for j = 1:nE
    for m = 1:d+1
        for n = 1:d+1 
            val = 0;
            for i = 1:d
                val = val+Vol_T(j)*Grads_T((d+1)*(j-1)+m,i)...
                    *u(d*(n4e(j,:)-1)+i)'*m_loc(:,n);
            end
            ctr = ctr+1;
            I(ctr) = n4e(j,m); J(ctr) = n4e(j,n); X_W(ctr) = val;
        end
    end
end
w = sparse(I,J,X_W,nC,nC); 
W = sparse(d*nC,d*nC);
for p = 1:d
    W(p:d:end,p:d:end) = w;
end