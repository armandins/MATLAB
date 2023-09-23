function show_nedelec(c4n,n4e,element2edges,edgeInd,u)
d = size(c4n,2); nE = size(n4e,1); idx = [1,3,6]; 
u_T = zeros(nE,d); mp_T = zeros(nE,d);
for j = 1:nE
    X_T = [ones(1,d+1);c4n(n4e(j,:),:)'];
    grads_T = X_T\[zeros(1,d);eye(d)];
    mp_T(j,:) = sum(c4n(n4e(j,:),:),1)/(d+1);
    for m = 1:idx(d)
        m1 = edgeInd(m,1); m2 = edgeInd(m,2);
        s_m = 2*(n4e(j,m1)<n4e(j,m2))-1;
        vec = s_m*(grads_T(m2,:)-grads_T(m1,:))/(d+1);
        u_T(j,:) = u_T(j,:)+u(element2edges(j,m))*vec;
    end
end
if d == 2
    triplot(n4e,c4n(:,1),c4n(:,2),'k'); hold on;
    quiver(mp_T(:,1),mp_T(:,2),u_T(:,1),u_T(:,2)); hold off;
else
    quiver3(mp_T(:,1),mp_T(:,2),mp_T(:,3),...
        u_T(:,1),u_T(:,2),u_T(:,3)); 
end

