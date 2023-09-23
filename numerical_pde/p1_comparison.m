function p1_comparison(d,red,assembly,solver)
[c4n,n4e,Db,Nb] = triang_cube(d); Db = [Db;Nb]; Nb = [];
for j = 1:red
    [c4n,n4e,Db,Nb] = red_refine(c4n,n4e,Db,Nb);
end
nC = size(c4n,1); h = 1/nC^(1/d);
dNodes = unique(Db); fNodes = setdiff(1:nC,dNodes);
[Vol_T,Grads_T,Mp_T] = nodal_basis(c4n,n4e);
Z = (1/(d+1))*Vol_T.*f(Mp_T); ZZ = repmat(Z,1,d+1);
switch assembly
    case 0
        [I,J,X] = fe_matrix_loop(c4n,n4e,Vol_T,Grads_T);
    case 1
        [I,J,X] = fe_matrix_vectorized_1(c4n,n4e,Vol_T,Grads_T);
    case 2
        [I,J,X] = fe_matrix_vectorized_2(c4n,n4e,Vol_T,Grads_T);
    case 3
        [I,J,X] = fe_matrix_mex(c4n,n4e,Vol_T,Grads_T);
end
s = sparse(I,J,X,nC,nC); u = u_D(c4n);
b = accumarray(n4e(:),ZZ(:),[nC,1])-s*u;
s_fN = s(fNodes,fNodes); b_fN = b(fNodes);
switch solver
    case 0
        u(fNodes) = s_fN\b_fN;
    case 1
        K = size(fNodes,2); eps_stop = h;
        C = spdiags(diag(s_fN),0,K,K);
        u(fNodes) = pcg(s_fN,b_fN,eps_stop,K,C);
end
show_p1(c4n,n4e,Db,Nb,u);
 
function val = f(x); val = ones(size(x,1),1);
function val = u_D(x); val = zeros(size(x,1),1);

function [Vol_T,Grads_T,Mp_T] = nodal_basis(c4n,n4e)
d = size(c4n,2); nE = size(n4e,1);
Grads_T = zeros((d+1)*nE,d); 
Vol_T = zeros(nE,1); Mp_T = zeros(nE,d);
for j = 1:nE
    X_T = [ones(1,d+1);c4n(n4e(j,:),:)'];
    Grads_T((j-1)*(d+1)+(1:d+1),:) = X_T\[zeros(1,d);eye(d)];
    Vol_T(j) = det(X_T)/factorial(d);
    Mp_T(j,:) = sum(c4n(n4e(j,:),:),1)/(d+1);
end