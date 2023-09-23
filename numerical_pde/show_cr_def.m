function show_cr_def(c4n,n4e,s4e,u)
nE = size(n4e,1); d = size(c4n,2);
Signum = ones(d+1)-d*eye(d+1); X = zeros((d+1)*nE,d); 
E = reshape(1:3*nE,3,nE)'; n4e_t = n4e';
for k = 1:d
    u_p1_disc_k = Signum*u(d*(s4e-1)+k)';
    X(:,k) = c4n(n4e_t(:),k)+u_p1_disc_k(:);
end
if d == 2
    trisurf(E,X(:,1),X(:,2),zeros((d+1)*nE,1)); view(0,90);
else
    for j = 1:nE
        tetramesh(1:d+1,X((d+1)*(j-1)+(1:d+1),:),0); hold on;
    end
    hold off; s = .1; axis([-s 1+s -s 1+s -s 1+s]); view(44,14)
end
