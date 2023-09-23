function multigrid(d_tmp,L)
global P d; d = d_tmp;
[c4n,n4e,Db,Nb] = triang_cube(d); Db = [Db;Nb]; Nb = [];
nC = size(c4n,1); fNodes_prev = setdiff(1:nC,unique(Db));
for ell = 1:L
    [c4n,n4e,Db,Nb,~,P1] = red_refine(c4n,n4e,Db,Nb);
    nC = size(c4n,1); fNodes = setdiff(1:nC,unique(Db));
    P{ell} = P1(fNodes,fNodes_prev);
    fNodes_prev = fNodes;
end
u = zeros(nC,1);
[s,m] = fe_matrices(c4n,n4e);
b = m(fNodes,:)*f(c4n);
A = s(fNodes,fNodes);
u(fNodes) = MG(A,b,L);
show_p1(c4n,n4e,Db,Nb,u)

function u = MG(A,b,ell)
global P;
nu_pre = 2; nu_post = 2; ell_0 = 2;
if ell == ell_0
    u = A\b;
else
    u_ini = zeros(size(b,1),1);
    u = richardson(A,b,u_ini,nu_pre,ell);
    A_coarse = P{ell}'*A*P{ell};
    r_coarse = P{ell}'*(b-A*u);
    c = MG(A_coarse,r_coarse,ell-1);
    u = u+P{ell}*c;
    u = richardson(A,b,u,nu_post,ell);
end

function u = richardson(A,b,u,nu,ell)
global d;
h = 2^(-ell); omega = h^(d-2)/10; 
for k = 1:nu
    u = u-omega*(A*u-b);
end

function val = f(x)
val = ones(size(x,1),1);