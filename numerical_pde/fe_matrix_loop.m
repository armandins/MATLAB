function [I,J,X] = fe_matrix_loop(c4n,n4e,Vol_T,Grads_T)
d = size(c4n,2); nE = size(n4e,1);
ctr = 0; ctr_max = (d+1)^2*nE;
I = zeros(ctr_max,1); J = zeros(ctr_max,1); 
X = zeros(ctr_max,1); 
for j = 1:nE
    grads_T = Grads_T((j-1)*(d+1)+(1:d+1),:);
    vol_T = Vol_T(j);
    for m = 1:d+1
        for n = 1:d+1
            ctr = ctr+1; I(ctr) = n4e(j,m); J(ctr) = n4e(j,n);
            X(ctr) = vol_T*grads_T(m,:)*grads_T(n,:)';
        end
    end
end