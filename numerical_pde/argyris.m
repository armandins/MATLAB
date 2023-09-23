function argyris(red)
theta = 1/2;
[c4n,n4e,Db,Nb] = triang_cube(2); Db = [Db;Nb];
for j = 1:red
    [c4n,n4e,Db,Nb] = red_refine(c4n,n4e,Db,Nb);
end
[s4e,sign_s4e,n4s,s4Db] = sides(n4e,Db,Nb);
nE = size(n4e,1); nS = size(n4s,1); nC = size(c4n,1);
dNodes_tmp = 6*(reshape(unique(Db),length(unique(Db)),1)-1);
dNodes = repmat(dNodes_tmp,1,3)+ones(size(dNodes_tmp,1),1)*(1:3);
fNodes = setdiff(1:6*nC+nS,[dNodes(:);6*nC+s4Db]);
A = sparse(6*nC+nS,6*nC+nS);
u = zeros(6*nC+nS,1); b = zeros(6*nC+nS,1);
for j = 1:nE
    loc_c4n = c4n(n4e(j,:),:);
    N = argyris_monomials_at_dofs(loc_c4n);    % conditioning ok?
    E = eye(21); E(19:21,19:21) = diag(sign_s4e(j,:));
    C = N\E;
    [xi,kappa] = argyris_gauss_quad_ref_deg_6(loc_c4n);
    X_T = [ones(1,3);loc_c4n'];
    [p,p_xx,p_yy,p_xy] = argyris_monomials_at_qps(xi);
    f_qp = f(xi);
    A_loc = zeros(21,21);
    b_loc = zeros(21,1);
    for i = 1:size(kappa,1);
        A_loc = A_loc+kappa(i)*det(X_T)*...
            ((p_xx(i,:)'*p_xx(i,:)+p_yy(i,:)'*p_yy(i,:)...
            +p_xx(i,:)'*p_yy(i,:)+p_yy(i,:)'*p_xx(i,:))...
            +(1-theta)*(2*p_xy(i,:)'*p_xy(i,:)...
            -p_xx(i,:)'*p_yy(i,:)-p_xx(i,:)'*p_yy(i,:)));
        b_loc = b_loc+kappa(i)*det(X_T)*f_qp(i)*p(i,:)';            
    end
    I_tmp = 6*(n4e(j,:)-1)';
    I = repmat(I_tmp,1,6)+ones(3,1)*(1:6);
    I = [I(:);6*nC+s4e(j,:)'];
    A(I,I) = A(I,I)+C'*A_loc*C;
    b(I) = b(I)+C'*b_loc;
end
u(fNodes) = A(fNodes,fNodes)\b(fNodes);
u_p1 = u(1:6:6*nC);
trisurf(n4e,c4n(:,1),c4n(:,2),u_p1);

function val = f(x)
val = 10*ones(size(x,1),1); 


