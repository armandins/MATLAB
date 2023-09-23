function streamline_diffusion(red)
d = 2; eps = 1e-5;
[c4n,n4e,Nb,Db] = triang_cube(d);
for j = 1:red
    [c4n,n4e,Db,Nb] = red_refine(c4n,n4e,Db,Nb);
end
[nC,d] = size(c4n); nE = size(n4e,1);
dNodes = unique(Db); fNodes = setdiff(1:nC,dNodes);
u = zeros(nC,1); tu_D = zeros(nC,1); 
ctr = 0; ctr_max = (d+1)^2*nE;
I = zeros(ctr_max,1); J = zeros(ctr_max,1); X = zeros(ctr_max,1);
for j = 1:nE
    X_T = [ones(1,d+1);c4n(n4e(j,:),:)'];
    grads_T = X_T\[zeros(1,d);eye(d)];
    vol_T = det(X_T)/factorial(d);
    mp_T = sum(c4n(n4e(j,:),:),1)/(d+1);
    b_T = b_field(mp_T);
    delta_T = vol_T^(1/d);
    for m = 1:d+1
        for n = 1:d+1
           ctr = ctr+1; I(ctr) = n4e(j,m); J(ctr) = n4e(j,n);
           X(ctr) = vol_T*(eps*grads_T(m,:)*grads_T(n,:)'...
               +delta_T*(b_T*grads_T(m,:)')*(b_T*grads_T(n,:)')...
               +b_T*grads_T(n,:)'/(d+1));
        end
    end
end
s = sparse(I,J,X,nC,nC);
for j = 1:nC
    tu_D(j) = u_D(c4n(j,:));
end
b = -s*tu_D; u(fNodes) = s(fNodes,fNodes)\b(fNodes); u = u+tu_D;
show_p1(c4n,n4e,Db,Nb,u);

function val = b_field(x)
[phi,~] = cart2pol(x(1),x(2)); 
val = [sin(phi),-cos(phi)];

function val = u_D(x)
val = 0; 
if (x(1)==0 && x(2)<=1/2); 
    val = 1; 
end