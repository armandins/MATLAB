function [I,J,X] = fe_matrix_vectorized_1(c4n,n4e,Vol_T,Grads_T)
d = size(c4n,2); nE = size(n4e,1); ctr_max = (d+1)^2*nE;
I = zeros(ctr_max,1); J = zeros(ctr_max,1); X = zeros(ctr_max,1);
for m = 1:d+1
    for n = 1:d+1
        idx = ((m-1)*(d+1)+(n-1))*nE+(1:nE);
        vals = Vol_T.*...
            sum(Grads_T(m:d+1:end,:).*Grads_T(n:d+1:end,:),2);
        I(idx) = n4e(:,m); J(idx) = n4e(:,n); X(idx) = vals;
    end
end
