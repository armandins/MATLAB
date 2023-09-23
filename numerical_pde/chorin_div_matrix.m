function [D,Grads_T,Vol_T] = chorin_div_matrix(c4n,n4e)
[nC,d] = size(c4n); nE = size(n4e,1);
ctr = 0; ctr_max = d*(d+1)^2*nE;
I = zeros(ctr_max,1); J = zeros(ctr_max,1); 
X_D = zeros(ctr_max,1); 
Vol_T =zeros(nE,1);
Grads_T = zeros((d+1)*nE,d);
for j = 1:nE
    X_T = [ones(1,d+1);c4n(n4e(j,:),:)'];
    grads_T = X_T\[zeros(1,d);eye(d)];
    vol_T = det(X_T)/factorial(d);
    for m = 1:d+1
        for n = 1:d+1
            for p = 1:d
                ctr = ctr+1; 
                I(ctr) = n4e(j,m); J(ctr) = d*(n4e(j,n)-1)+p;
                X_D(ctr) = vol_T*grads_T(m,p)'/(d+1);
            end
        end
    end
    Vol_T(j) = vol_T;
    Grads_T((d+1)*(j-1)+(1:d+1),:) = grads_T;
end
D = sparse(I,J,X_D,nC,d*nC); 