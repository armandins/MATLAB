function gen_triang_cyl_w_hole(d_tmp)
addpath('~/auxiliary/distmesh');
global d r_sph L_cyl r_cyl;
d = d_tmp; L_cyl = 2; r_cyl = 1; r_sph = 1/2;
R = 2; h_min = 0.1; fixed = [];
box = [-R*ones(1,d);R*ones(1,d)];
[c4n,n4e] = distmeshnd(@s,@phi,h_min,box,fixed);
n4e = fix_orientation(c4n,n4e);
Db = find_bdy_sides(n4e); Nb = [];
str = strcat('save triang_cyl_w_hole_',num2str(d),...
    'd.mat c4n n4e Db Nb');
% eval(str);

function val = s(x)
global d r_sph L_cyl r_cyl;
dist_hor = abs(x(:,1))-L_cyl;
dist_rad = sqrt(sum(x(:,2:d).^2,2))-r_cyl;
dist_cyl = max(dist_hor,dist_rad);
idx = dist_hor>0 & dist_rad>0;
dist_cyl(idx) = sqrt(dist_hor(idx).^2+dist_rad(idx).^2);
dist_compl_sph = r_sph-sqrt(sum(x.^2,2));
val = max(dist_cyl,dist_compl_sph);

function val = phi(x)
global r_sph;
dist_sph = sqrt(sum(x.^2,2))-r_sph;
val = min(dist_sph+1,2);

function bdy = find_bdy_sides(n4e)
d = size(n4e,2)-1;
if d == 2
    all_sides = [n4e(:,[1,2]);n4e(:,[2,3]);n4e(:,[3,1])];
elseif d == 3
    all_sides = [n4e(:,[2,4,3]);n4e(:,[1,3,4]);...
        n4e(:,[1,4,2]);n4e(:,[1,2,3])];
end
[sides,~,j] = unique(sort(all_sides,2),'rows');
valence = accumarray(j(:),1);
bdy = sides(valence==1,:);

function n4e = fix_orientation(c4n,n4e)
global d; nE = size(n4e,1); or_Vol_T = zeros(nE,1); 
for j = 1:nE
    X_T = [ones(1,d+1);c4n(n4e(j,:),:)'];
    or_Vol_T(j) = det(X_T)/factorial(d);
end
n4e(or_Vol_T<0,[1,2]) = n4e(or_Vol_T<0,[2,1]);
