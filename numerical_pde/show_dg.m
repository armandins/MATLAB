function show_dg(c4n,n4e,u)
nE = size(n4e,1); d = size(c4n,2);
u_elements = reshape(u,d+1,nE)';
u_mp = sum(u_elements,2)/(d+1);
Signum = ones(d+1)-d*eye(d+1);
if d == 2
    E = reshape(1:3*nE,3,nE)'; n4e_t = n4e';
    X = c4n(n4e_t(:),1); Y = c4n(n4e_t(:),2); Z = Signum*u_elements';
    trisurf(E,X,Y,Z(:)); 
else
    tetramesh(n4e,c4n,u_mp);
end