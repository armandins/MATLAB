function show_p2_iso(c4n,n4e,u)
nE = size(n4e,1); fine = 5;
c4n_ref = [0 0;1 0;0 1]; n4e_ref = [1 2 3]; 
Db_ref = [1 2;2 3;3 1]; 
for j = 1:fine
    [c4n_ref,n4e_ref,Db_ref] = ...
        red_refine(c4n_ref,n4e_ref,Db_ref,[]);
end
r = c4n_ref(:,1); s = c4n_ref(:,2);
psi = [1-r-s,r,s,4*r.*(1-r-s),4*r.*s,4*s.*(1-r-s)];
N = [1,1,0;0,1,1;1,0,1]/2;
hold off;
for j = 1:nE
    K_T = find(n4e(j,:));
    P = zeros(6,2);
    P(K_T,:) = c4n(n4e(j,K_T),:);
    P(4:6,:) = P(4:6,:)+((n4e(j,4:6)==0)'*[1,1]).*(N*P(1:3,:));
    D(1:3,:) = P(1:3,:);
    D(4:6,:) = P(4:6,:)-(N*P(1:3,:));
    c4n_def_T = psi*D;
    U = zeros(6,1);
    U(K_T) = u(n4e(j,K_T));
    u_def = psi*U;
    trisurf(n4e_ref,c4n_def_T(:,1),c4n_def_T(:,2),u_def); 
    shading flat; hold on;
    for k = 1:size(Db_ref,1);
        X = c4n_def_T(Db_ref(k,:),1);
        Y = c4n_def_T(Db_ref(k,:),2);
        Z = u_def(Db_ref(k,:));
        plot3(X,Y,Z,'k');
    end
end
view(22,36); hold off;
