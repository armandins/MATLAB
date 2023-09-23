function p1_poisson(d,red)
[c4n,n4e,Db,Nb] = triang_cube(d);
for j = 1:red
    [c4n,n4e,Db,Nb,P0,P1] = red_refine(c4n,n4e,Db,Nb);
end
[nC,d] = size(c4n); nE = size(n4e,1); nNb = size(Nb,1);
dNodes = unique(Db); fNodes = setdiff(1:nC,dNodes);
u = zeros(nC,1); tu_D = zeros(nC,1); b = zeros(nC,1); 
ctr = 0; ctr_max = (d+1)^2*nE;
I = zeros(ctr_max,1); J = zeros(ctr_max,1); X = zeros(ctr_max,1);
for j = 1:nE
    X_T = [ones(1,d+1);c4n(n4e(j,:),:)'];
    grads_T = X_T\[zeros(1,d);eye(d)];
    vol_T = det(X_T)/factorial(d);
    mp_T = sum(c4n(n4e(j,:),:),1)/(d+1);
    for m = 1:d+1
        b(n4e(j,m)) = b(n4e(j,m))+(1/(d+1))*vol_T*f(mp_T);
        for n = 1:d+1
            ctr = ctr+1; I(ctr) = n4e(j,m); J(ctr) = n4e(j,n);
            X(ctr) = vol_T*grads_T(m,:)*grads_T(n,:)';
        end
    end
end
s = sparse(I,J,X,nC,nC);
for j = 1:nNb
    if d == 1
        vol_S = 1;
    elseif d == 2
        vol_S = norm(c4n(Nb(j,1),:)-c4n(Nb(j,2),:));
    elseif d == 3
        vol_S = norm(cross(c4n(Nb(j,3),:)-c4n(Nb(j,1),:),...
            c4n(Nb(j,2),:)-c4n(Nb(j,1),:)),2)/2;
    end
    mp_S = sum(c4n(Nb(j,:),:),1)/d;
    for k = 1:d
        b(Nb(j,k)) = b(Nb(j,k))+(1/d)*vol_S*g(mp_S);
    end
end
for j = 1:nC
    tu_D(j) = u_D(c4n(j,:));
end
b = b-s*tu_D; u(fNodes) = s(fNodes,fNodes)\b(fNodes); u = u+tu_D;
if d == 1; plot(c4n(n4e),u(n4e));
elseif d == 2; trisurf(n4e,c4n(:,1),c4n(:,2),u);
elseif d == 3; trisurf([Db;Nb],c4n(:,1),c4n(:,2),c4n(:,3),u);
end

function val = f(x); val = 1;
function val = g(x); val = 1;
function val = u_D(x); val = sin(2*pi*x(:,1));
