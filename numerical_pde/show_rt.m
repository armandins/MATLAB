function show_rt(c4n,n4e,u,p,s4e,sign_s4e)
nE = size(n4e,1);
u_T = zeros(nE,2); mp_T = zeros(nE,2);
shift1 = [2,3,1]; shift2 = [3,1,2];
for j = 1:nE
    area_T = det([1,1,1;c4n(n4e(j,:),:)'])/2;
    mp_T(j,:) = sum(c4n(n4e(j,:),:))/3;
    for m = 1:3
        side_m = s4e(j,m);
        length_m = ...
            norm(c4n(n4e(j,shift1(m)),:)-c4n(n4e(j,shift2(m)),:));
        signum_m = sign_s4e(j,m);
        u_T(j,:) = u_T(j,:)+u(side_m)*signum_m*length_m...
            *(c4n(n4e(j,m),:)-mp_T(j,:))/(2*area_T);
    end
end
E = reshape(1:3*nE,3,nE)'; n4e_t = n4e'; 
X = c4n(n4e_t(:),1); Y = c4n(n4e_t(:),2); P = repmat(p,1,3)';
figure(1); clf; subplot(1,2,1); trisurf(E,X,Y,P(:));
subplot(1,2,2); triplot(n4e,c4n(:,1),c4n(:,2),':k'); hold on;
quiver(mp_T(:,1),mp_T(:,2),u_T(:,1),u_T(:,2));