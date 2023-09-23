function p1_graded(N,beta)
[c4n,n4e,Db] = Lshape_graded(N,beta);
nC = size(c4n,1); 
dNodes = unique(Db); fNodes = setdiff(1:nC,dNodes);
[s,m] = fe_matrices(c4n,n4e);
u = zeros(nC,1); b = m*f(c4n);
u(fNodes) = s(fNodes,fNodes)\b(fNodes);
show_p1(c4n,n4e,Db,[],u);

function [c4n,n4e,Db] = Lshape_graded(N,beta) 
c4n_macro = [-1 -1;0 -1;-1 0;0 0;1 0;-1 1;0 1;1 1];
n4e_macro = [4 5 8;4 8 7;4 7 6;4 6 3;4 3 1;4 1 2];
nE_macro = size(n4e_macro,1);
[c4n_micro,n4e_micro] = grad_grid_ref(N,beta);
nC_mi = size(c4n_micro,1);
n4e = []; c4n = [];
for j = 1:nE_macro
    phi_0 = (1-c4n_micro(:,1)-c4n_micro(:,2));
    phi_1 = c4n_micro(:,1); phi_2 = c4n_micro(:,2); 
    c4n_transf = phi_0*c4n_macro(n4e_macro(j,1),:)...
        +phi_1*c4n_macro(n4e_macro(j,2),:)...
        +phi_2*c4n_macro(n4e_macro(j,3),:);
    n4e = [n4e;n4e_micro+(j-1)*nC_mi];
    c4n = [c4n;c4n_transf];
end
[c4n,~,j] = unique(c4n,'rows');
n4e = j(n4e);
Db = find_bdy_sides(n4e);

function [c4n,n4e] = grad_grid_ref(N,beta)
c4n = [0,0]; n4e = [1 2 3];
for j = 1:N
    xi = (j/N)^beta;
    c4n = [c4n;ones(j+1,1)*[xi,0]+(0:j)'/j*[-xi,xi]];
end
for j = 1:N-1
    n4e = [n4e;j*(j+1)/2+[(1:j)',j+2+(1:j)',1+(1:j)']];
    n4e = [n4e;j*(j+1)/2+[(1:j+1)',j+1+(1:j+1)',j+2+(1:j+1)']];
end

function bdy = find_bdy_sides(n4e)
all_sides = [n4e(:,[1,2]);n4e(:,[2,3]);n4e(:,[3,1])];
[sides,~,j] = unique(sort(all_sides,2),'rows');
valence = accumarray(j(:),1);
bdy = sides(valence==1,:);

function val = f(x); val = ones(size(x,1),1);