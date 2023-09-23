function maxwell(d,red)
[c4n,n4e,Db,Nb] = triang_cube(d); o_sq = 1;
for j = 1:red
    [c4n,n4e,Db,Nb] = red_refine(c4n,n4e,Db,Nb);
end
nE = size(n4e,1); d = size(c4n,2); 
switch d
    case 2; edgeInd = [1 2,1 3,2 3];
    case 3; edgeInd = [1 2,1 3,1 4,2 3,2 4,3 4];
end
[edges,el2edges,Db2edges] = edge_data(n4e,Db,Nb);
nEdges = size(edges,1); fEdges = setdiff(1:nEdges,Db2edges);
u = zeros(nEdges,1); b = zeros(nEdges,1);
idx = [1,3,6]; edgeInd = reshape(edgeInd,2,idx(d))';
m_loc = (eye(d+1)+ones(d+1,d+1))/((d+2)*(d+1));
max_ctr = idx(d)*nE; ctr = 0;
I = zeros(max_ctr,1); J = zeros(max_ctr,1);
X_C = zeros(max_ctr,1); X_M = zeros(max_ctr,1);
for j = 1:nE
    X_T = [ones(1,d+1);c4n(n4e(j,:),:)'];
    grads_T = X_T\[zeros(1,d);eye(d)];
    vol_T = det(X_T)/factorial(d);
    mp_T = sum(c4n(n4e(j,:),:),1)/(d+1);
    for m = 1:idx(d)
       m1 = edgeInd(m,1); m2 = edgeInd(m,2);
       s_m = 2*(n4e(j,m1)<n4e(j,m2))-1;
       curl_psi_m = 2*s_m*wedge(grads_T(m1,:),grads_T(m2,:),d);
       for n = 1:idx(d)
          ctr = ctr+1;
          I(ctr) = el2edges(j,m); J(ctr) = el2edges(j,n);
          n1 = edgeInd(n,1); n2 = edgeInd(n,2);
          s_n = 2*(n4e(j,n1)<n4e(j,n2))-1;
          curl_psi_n = 2*s_n*wedge(grads_T(n1,:),grads_T(n2,:),d);
          X_C(ctr) = vol_T*dot(curl_psi_m,curl_psi_n);
          X_M(ctr) = vol_T*s_m*s_n*...
              (m_loc(m1,n1)*grads_T(m2,:)*grads_T(n2,:)'...
              -m_loc(m2,n1)*grads_T(m1,:)*grads_T(n2,:)'...
              -m_loc(m1,n2)*grads_T(m2,:)*grads_T(n1,:)'...
              +m_loc(m2,n2)*grads_T(m1,:)*grads_T(n1,:)');
       end
       b(el2edges(j,m))= b(el2edges(j,m))+vol_T*f(mp_T,o_sq,d)'...
           *s_m*(grads_T(m2,:)-grads_T(m1,:))'/(d+1);
    end
end
C = sparse(I,J,X_C); M = sparse(I,J,X_M); X = C-o_sq*M;
u(fEdges) = X(fEdges,fEdges)\b(fEdges);
show_nedelec(c4n,n4e,el2edges,edgeInd,u);

function val = f(x,o_sq,d)
val = (pi^2-o_sq)*[sin(pi*x(2));sin(pi*x(1));zeros(d-2,1)];

function val = wedge(v,w,d)
if d == 2; val = det([v;w]); else val = cross(v,w); end
