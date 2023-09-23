function [I,J,X] = fe_matrix_vectorized_2(c4n,n4e,Vol_T,Grads_T)
d = size(c4n,2); nE = size(n4e,1);
I = repmat(reshape(n4e',(d+1)*nE,1),1,d+1)';
J = repmat(n4e,1,d+1)';
B = reshape(repmat(Grads_T',d+1,1),d,(d+1)^2*nE)';
C = reshape(repmat(reshape(Grads_T',(d+1)*d,nE),d+1,1),...
    d,(d+1)^2*nE)';
rep_Vol_T = repmat(Vol_T,1,(d+1)^2)';
X = rep_Vol_T(:).*sum(B.*C,2);
I = I(:); J = J(:); X = X(:);
