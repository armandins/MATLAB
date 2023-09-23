function stokes_cr
d = 2; mu = 1; 
load triang_cyl_w_hole_2d;
idx = (min([c4n(Db(:,1),1),c4n(Db(:,2),1)],[],2)>=2.0);
Nb = Db(idx,:); Db(idx,:) = [];
[s4e,~,n4s,s4Db,s4Nb] = sides(n4e,Db,Nb);
nS = size(n4s,1); nE = size(n4e,1); nNb = size(Nb,1); 
fNodes = setdiff(1:2*nS+nE,[d*(s4Db-1)+1;d*(s4Db-1)+2]);
u = zeros(d*nS,1); tu_D = u; mp_T = zeros(nE,d);
x = zeros(d*nS+nE,1); tx = x; b = x;
ctr_A = 0; ctr_A_max = (d+1)^2*nE;
X = zeros(ctr_A_max,1); I = X; J = X;
ctr_B = 0; ctr_B_max = (d+1)*nE;
Y = zeros(ctr_B_max,1); K = Y; L = Y;
for j = 1:nE
    X_T = [ones(1,d+1);c4n(n4e(j,:),:)'];
    grads_T = X_T\[zeros(1,d);eye(d)];
    vol_T = det(X_T)/factorial(d);
    mp_T(j,:) = sum(c4n(n4e(j,:),:),1)/(d+1);
    for m = 1:d+1
      b(d*(s4e(j,m)-1)+(1:d)) = b(d*s4e(j,m)-1)...
          +(1/(d+1))*vol_T*f(mp_T(j,:));
      for n = 1:d+1
        ctr_A = ctr_A+1; I(ctr_A) = s4e(j,m); J(ctr_A) = s4e(j,n);
        X(ctr_A) = mu*d^2*vol_T*grads_T(m,:)*grads_T(n,:)';
      end
      K(ctr_B+(1:d)) = j; L(ctr_B+(1:d)) = d*(s4e(j,m)-1)+(1:d);
      Y(ctr_B+(1:d)) = d*vol_T*grads_T(m,:);
      ctr_B = ctr_B+d;
    end
end
A = sparse([2*I-1;2*I],[2*J-1;2*J],[X;X],d*nS,d*nS);      % d=2
B = sparse(K,L,Y,nE,d*nS);
for j = 1:nNb
    vol_S = norm(c4n(Nb(j,1),:)-c4n(Nb(j,2),:));          % d=2
    mp_S = sum(c4n(n4s(j,:),:),1)/d;
    b(d*(s4Nb(j)-1)+(1:d)) = b(d*(s4Nb(j)-1)+(1:d))+vol_S*g(mp_S);
end
for j = 1:nS
    mp_S = sum(c4n(n4s(j,:),:),1)/d;
    tu_D(d*(j-1)+(1:d)) = u_D(mp_S);
end
tx(1:d*nS) = tu_D;
G = [A,B';B,sparse(nE,nE)]; b = b-G*tx;
x(fNodes) = G(fNodes,fNodes)\b(fNodes);
u = x(1:d*nS); p = x(d*nS+(1:nE));
show_stokes_cr(c4n,n4e,s4e,u,p,mp_T)

function val = f(x); d = size(x,2); val = zeros(d,1);
function val = g(x); d = size(x,2); val = zeros(d,1);
function val = u_D(x); d = size(x,2); val = zeros(d,1);
if x(1)<=-2; val(1) = 1; end