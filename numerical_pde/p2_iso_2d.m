function p2_iso_2d
[c4n,n4e,Db] = triang_p2_iso();
nC = size(c4n,1); nE = size(n4e,1);
A = sparse(nC,nC); b = zeros(nC,1); u = zeros(nC,1); 
fNodes = setdiff(1:nC,unique(Db));
[phi,d1_phi,d2_phi,kappa] = quad_p2_iso();
N = [1,1,0;0,1,1;1,0,1]/2;
for j = 1:nE
    K_T = find(n4e(j,:));
    P = zeros(6,2);
    P(K_T,:) = c4n(n4e(j,K_T),:);
    P(4:6,:) = P(4:6,:)+((n4e(j,4:6)==0)'*[1,1]).*(N*P(1:3,:));
    D(1:3,:) = P(1:3,:);
    D(4:6,:) = P(4:6,:)-(N*P(1:3,:));
    M = zeros(6,6);
    det_D_Psi = zeros(1,size(kappa,2));
    for m = 1:size(kappa,2)
        D_Psi = [d1_phi(m,:);d2_phi(m,:)]*D;
        D_phi_transp = D_Psi\[d1_phi(m,:);d2_phi(m,:)];
        det_D_Psi(m) = abs(det(D_Psi));
        M = M+kappa(m)*(D_phi_transp'*D_phi_transp)*det_D_Psi(m);
    end
    A(n4e(j,K_T),n4e(j,K_T)) = ...
        A(n4e(j,K_T),n4e(j,K_T))+M(K_T,K_T);
    val_b = kappa.*det_D_Psi.*f(phi*D)'*phi;
    b(n4e(j,K_T)) = b(n4e(j,K_T))+val_b(K_T)';
end
idx = find(Db(:,3));
u(unique(Db(:,1:2))) = u_D(c4n(unique(Db(:,1:2)),:));
u(Db(idx,3))=u_D(c4n(Db(idx,3),:))-(u(Db(idx,1))+u(Db(idx,2)))/2;
b = b-A*u;
u(fNodes) = A(fNodes,fNodes)\b(fNodes);
show_p2_iso(c4n,n4e,u)

function val = f(x)
val = 0*ones(size(x,1),1);

function val = u_D(x)
[phi,r] = cart2pol(x(:,1),x(:,2));
phi = phi+2*pi*(phi<0);
val = r.^(2/3).*sin(2/3*phi);

function [phi,d1_phi,d2_phi,kappa] = quad_p2_iso()
pos = [6-sqrt(15),9+2*sqrt(15),6+sqrt(15),9-2*sqrt(15),7]/21;
r = pos([1,2,1,3,3,4,5])'; s = pos([1,1,2,4,3,3,5])';
wts = [155-sqrt(15),155+sqrt(15),270]/2400;
kappa = wts([1,1,1,2,2,2,3]);
one = ones(size(kappa,2),1);
phi = [1-r-s,r,s,4*r.*(1-r-s),4*r.*s,4*s.*(1-r-s)];
d1_phi = [-one,one,0*one,4*(1-2*r-s),4*s,-4*s];
d2_phi = [-one,0*one,one,-4*r,4*r,4*(1-r-2*s)];