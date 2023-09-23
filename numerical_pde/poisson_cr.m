function poisson_cr(d,red)
[c4n,n4e,Db,Nb] = triang_cube(d);
for j = 1:red
    [c4n,n4e,Db,Nb] = red_refine(c4n,n4e,Db,Nb);
end
[s4e,~,n4s,s4Db,s4Nb] = sides(n4e,Db,Nb);
nS = size(n4s,1); nE = size(n4e,1); nNb = size(Nb,1); 
fSides = setdiff(1:nS,s4Db);
u = zeros(nS,1); tu_D = zeros(nS,1); b = zeros(nS,1); 
ctr = 0; ctr_max = (d+1)^2*nE;
I = zeros(ctr_max,1); J = zeros(ctr_max,1); X = zeros(ctr_max,1);
for j = 1:nE
    X_T = [ones(1,d+1);c4n(n4e(j,:),:)'];
    grads_T = X_T\[zeros(1,d);eye(d)];
    vol_T = det(X_T)/factorial(d);
    mp_T = sum(c4n(n4e(j,:),:),1)/(d+1);
    for m = 1:d+1
        b(s4e(j,m)) = b(s4e(j,m))+(1/(d+1))*vol_T*f(mp_T);
        for n = 1:d+1
            ctr = ctr+1; I(ctr) = s4e(j,m); J(ctr) = s4e(j,n);
            X(ctr) = d^2*vol_T*grads_T(m,:)*grads_T(n,:)';
        end
    end
end
s = sparse(I,J,X,nS,nS);
for j = 1:nNb
    if d == 2
        vol_S = norm(c4n(Nb(j,1),:)-c4n(Nb(j,2),:));
    elseif d == 3
        vol_S = norm(cross(c4n(Nb(j,3),:)-c4n(Nb(j,1),:),...
            c4n(Nb(j,2),:)-c4n(Nb(j,1),:)),2)/2;
    end
    mp_S = sum(c4n(Nb(j,:),:),1)/d;
    for k = 1:d
        b(s4Nb(j)) = b(s4Nb(j))+(1/d)*vol_S*g(mp_S);
    end
end
for j = 1:nS
    mp_S = sum(c4n(n4s(j,:),:),1)/d;
    tu_D(j) = u_D(mp_S);
end
b = b-s*tu_D; u(fSides) = s(fSides,fSides)\b(fSides); u = u+tu_D;
show_cr(c4n,n4e,s4e,u)

function val = f(x); val = 1;
function val = g(x); val = 1;
function val = u_D(x); val = sin(2*pi*x(:,1));