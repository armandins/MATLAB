function chorin_projection(d_tmp,red)
global d; d = d_tmp;
tau = 2^(-red)/10; nu = .01; T = 10; K = ceil(T/tau);
str = strcat('load triang_cyl_w_hole_',num2str(d),'d');
eval(str); 
for j = 1:red
    [c4n,n4e,Db,Nb] = red_refine(c4n,n4e,Db,Nb);
end
nC = size(c4n,1); 
dNodes = unique(Db); fNodes = setdiff(1:nC,dNodes);
FNodes = repmat(d*(fNodes-1),d,1)+(1:d)'*ones(1,size(fNodes,2),1);
FNodes = FNodes(:);
[s,m] = fe_matrices(c4n,n4e);
S = sparse(d*nC,d*nC); M = sparse(d*nC,d*nC);
for j = 1:d
    S(j:d:d*nC,j:d:d*nC) = s; M(j:d:d*nC,j:d:d*nC) = m;
end
[D,Grads_T,Vol_T] = chorin_div_matrix(c4n,n4e);
u_old = u_D(0,c4n);
for k = 1:K
    t = k*tau
    tu_new = u_D(t,c4n); tu_new(FNodes) = 0;
    % W = chorin_conv_matrix(c4n,n4e,Grads_T,Vol_T,u_old);
    W = chorin_conv_matrix_vec(c4n,n4e,Grads_T,Vol_T,u_old);
    A = M+tau*nu*S+tau*W;
    b = tau*M*f(t,c4n)+M*u_old-A*tu_new;
    tu_new(FNodes) = A(FNodes,FNodes)\b(FNodes); 
    c = (1/tau)*D*(tu_new-u_D(t,c4n))-(1/tau)*m*div_u_D(t,c4n);
    p = zeros(nC,1); 
    p(2:nC) = s(2:nC,2:nC)\c(2:nC);
    Pi_nabla_p = M\(D'*p);
    u_new = tu_new-tau*Pi_nabla_p;
    show_chorin(c4n,n4e,u_new,p);
    u_old = u_new;
end

function val = u_D(t,x) 
global d; val = zeros(d,size(x,1)); idx = find(abs(x(:,1))>1);
val(1,idx) = sin(t)*(abs(x(idx,1))-1).*(sum(x(idx,2:d).^2,2)-1);
val = val(:);

function val = div_u_D(t,x)
global d; val = zeros(size(x,1),1); idx = find(abs(x(:,1))>1);
val(idx) = sin(t)*sign(x(idx,1)).*(sum(x(idx,2:d).^2,2)-1); 

function val = f(t,x)
global d; val = zeros(d*size(x,1),1);
