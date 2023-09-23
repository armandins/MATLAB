function W = chorin_conv_matrix_vec(c4n,n4e,Grads_T,Vol_T,u)
[nC,d] = size(c4n); nE = size(n4e,1); n4e_t = n4e';
m_loc = (ones(d+1,d+1)+eye(d+1))/((d+1)*(d+2));
X_w = zeros((d+1)^2,nE);
for i = 1:d
    X1 = repmat(Vol_T,1,(d+1)^2);
    X2 = reshape(repmat(Grads_T(:,i),1,d+1)',(d+1)^2,nE)';
    X3 = repmat(u(d*(n4e-1)+i)*m_loc,1,d+1);
    X_w = X_w+(X1.*X2.*X3)';
end
I = repmat(reshape(n4e_t,1,(d+1)*nE),d+1,1);
J = repmat(n4e_t,d+1,1);
if d == 2   
    W = sparse([2*I-1,2*I],[2*J-1,2*J],[X_w,X_w],d*nC,d*nC);
elseif d == 3
    W = sparse([3*I-2,3*I-1,3*I],[3*J-2,3*J-1,3*J],...
        [X_w,X_w,X_w],d*nC,d*nC);
end
