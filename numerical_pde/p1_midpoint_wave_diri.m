function p1_midpoint_wave_diri(d,red)
T = 10; alpha = 1;
[c4n,n4e,Db,Nb] = triang_cube(d); Db = [Db;Nb]; Nb = [];
for j = 1:red
    [c4n,n4e,Db,Nb] = red_refine(c4n,n4e,Db,Nb);
end
nC = size(c4n,1); h = 2^(-red); tau = h^alpha/4; K = floor(T/tau);
dNodes = unique(Db); fNodes = setdiff(1:nC,dNodes);
u_old_1 = u_0(c4n)+tau*v_0(c4n)+(tau^2/2)*(Del_u_0(c4n)+f(0,c4n));
u_old_2 = u_0(c4n); u_new = zeros(nC,1);
[s,m,m_lumped] = fe_matrices(c4n,n4e); 
for k = 2:K
    b = (1/tau^2)*m*(2*u_old_1-u_old_2)...
        -(1/4)*s*(2*u_old_1+u_old_2)+m*f((k-1)*tau,c4n);
    X = (1/tau^2)*m+(1/4)*s;    
    u_new(fNodes) = X(fNodes,fNodes)\b(fNodes);
    show_p1(c4n,n4e,Db,Nb,u_new); 
    axis([0,1,0,1,-1,1,-1,1]); pause(.05);
    u_old_2 = u_old_1; u_old_1 = u_new;
end

function val = f(t,x); val = -ones(size(x,1),1);
function val = u_0(x); val = zeros(size(x,1),1); 
function val = v_0(x); val = 4*sin(pi*x(:,1)).*sin(pi*x(:,2));
function val = Del_u_0(x); val = zeros(size(x,1),1);