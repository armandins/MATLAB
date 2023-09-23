function [Vol_T,Grads_T,Mp_T] = nodal_basis(c4n,n4e)
d = size(c4n,2); nE = size(n4e,1);
Grads_T = zeros((d+1)*nE,d); Vol_T = zeros(nE,1); Mp_T = zeros(nE,d);
for j = 1:nE
    X_T = [ones(1,d+1);c4n(n4e(j,:),:)'];
    Grads_T((j-1)*(d+1)+(1:d+1),:) = X_T\[zeros(1,d);eye(d)];
    Vol_T(j) = det(X_T)/factorial(d);
    Mp_T(j,:) = sum(c4n(n4e(j,:),:),1)/(d+1);
end
