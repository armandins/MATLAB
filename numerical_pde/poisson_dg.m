function poisson_dg(d,red)
global beta gamma;
sigma = 1; beta = 10; gamma = 1;
[c4n,n4e,Db,Nb] = triang_cube(d); Db = [Db;Nb]; Nb = [];
for j = 1:red
    [c4n,n4e,Db,Nb] = red_refine(c4n,n4e,Db,Nb);
end
[s4e,~,n4s,~,~,e4s] = sides(n4e,Db,Nb);
nS = size(n4s,1); nE = size(n4e,1); 
b = zeros((d+1)*nE,1); 
ctr = 0; ctr_max = (d+1)^2*nE;
grads = zeros((d+1)*nE,d); 
normals_S = zeros(nS,d); vol_S = zeros(nS,1);
I = zeros(ctr_max,1); J = zeros(ctr_max,1); X = zeros(ctr_max,1);
for j = 1:nE
    X_T = [ones(1,d+1);c4n(n4e(j,:),:)'];
    grads_T = X_T\[zeros(1,d);eye(d)];
    vol_T = det(X_T)/factorial(d);
    mp_T = sum(c4n(n4e(j,:),:),1)/(d+1);
    for m = 1:d+1
        b((j-1)*(d+1)+m) = (1/(d+1))*vol_T*f(mp_T);
        for n = 1:d+1
            ctr = ctr+1; 
            I(ctr) = (j-1)*(d+1)+m; J(ctr) = (j-1)*(d+1)+n;
            X(ctr) = d^2*vol_T*grads_T(m,:)*grads_T(n,:)';
        end
    end
    grads((j-1)*(d+1)+(1:d+1),:) = grads_T;
    heights = 1./sqrt(sum(grads_T.^2,2));
    vol_S(s4e(j,:)) = factorial(d)*vol_T./heights;
    normals_S(s4e(j,:),:) = -grads_T.*(heights*ones(1,d)); 
end
s_elements = sparse(I,J,X,(d+1)*nE,(d+1)*nE);
[s_sides,s_penal] = ...
    dg_side_matrices(n4e,e4s,n4s,vol_S,grads,normals_S);
s = s_elements+s_sides+sigma*s_sides'+s_penal;
u = s\b; 
show_dg(c4n,n4e,u)

function val = f(x); val = 1;
