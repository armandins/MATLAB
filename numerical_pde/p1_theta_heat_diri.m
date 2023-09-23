function p1_theta_heat_diri(d,red)
T = 10; theta = 1/2; alpha = 1;
[c4n,n4e,Db,Nb] = triang_cube(d); Db = [Db;Nb]; Nb = [];
for j = 1:red
    [c4n,n4e,Db,Nb] = red_refine(c4n,n4e,Db,Nb);
end
nC = size(c4n,1); h = 2^(-red); tau = h^alpha/4; K = floor(T/tau);
dNodes = unique(Db); fNodes = setdiff(1:nC,dNodes);
u_old = u_0(c4n); u_new = zeros(nC,1);
[s,m,m_lumped] = fe_matrices(c4n,n4e); 
for k = 1:K
    t_k_theta = (k-1+theta)*tau;
    b = (1/tau)*m*u_old-(1-theta)*s*u_old+m*f(t_k_theta,c4n);
    X = (1/tau)*m+theta*s;    
    u_new(fNodes) = X(fNodes,fNodes)\b(fNodes);
    show_p1(c4n,n4e,Db,Nb,u_new);
    axis([0,1,0,1,0,.25,0,.25]); pause(.1);
    u_old = u_new;
end

function val = f(t,x); val = ones(size(x,1),1);
function val = u_0(x); val = sin(2*pi*x(:,1)).*sin(2*pi*x(:,2)); 
