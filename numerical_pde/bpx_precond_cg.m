function bpx_precond_cg(d_tmp,L_tmp)
global h d L P_full; d = d_tmp; L = L_tmp;
[c4n,n4e,Db,Nb] = triang_cube(d); Db = [Db;Nb]; Nb = [];
nC = size(c4n,1); fNodes_prev = setdiff(1:nC,unique(Db));
h = zeros(L,1);
for ell = 1:L
    [c4n,n4e,Db,Nb,~,P1] = red_refine(c4n,n4e,Db,Nb);
    nC = size(c4n,1); fNodes = setdiff(1:nC,unique(Db));
    P{ell} = P1(fNodes,fNodes_prev);
    fNodes_prev = fNodes; h(ell) = 2^(-ell);
end
nfNodes = size(fNodes,2);
P_full{L} = speye(nfNodes);
for ell = L-1:-1:1
    P_full{ell} = P_full{ell+1}*P{ell+1};
end
[s,m] = fe_matrices(c4n,n4e);
A = s(fNodes,fNodes); 
b = m(fNodes,:)*f(c4n); 
u = zeros(nC,1);
u(fNodes) = cg_precond(u(fNodes),A,b);
show_p1(c4n,n4e,Db,Nb,u)

function x = cg_precond(x,A,b)
r = b-A*x; z = apply_bpx(r); d = z; 
rz_old = r'*z; eps = 1e-4;
while sqrt(r'*r) > eps
    alpha = rz_old/(d'*A*d);
    x = x+alpha*d; 
    r = r-alpha*A*d;
    z = apply_bpx(r); 
    rz_new = r'*z;
    beta = rz_new/rz_old;
    d = z+beta*d; 
    rz_old = rz_new;
end

function Cr = apply_bpx(r)
global h d L P_full;
Cr = zeros(size(r));
for ell = 1:L
    Cr = Cr+h(ell)^(2-d)*P_full{ell}*(P_full{ell}'*r);
end

function val = f(x); val = ones(size(x,1),1);