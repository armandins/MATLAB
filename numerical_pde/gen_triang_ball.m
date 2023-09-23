function [c4n,n4e] = gen_triang_ball(d_tmp,r_tmp,h_min,alpha_tmp)
addpath('~/auxiliary/distmesh');
global d r alpha;
d = d_tmp; r = r_tmp; alpha = alpha_tmp;
R = r; fixed = [];
box = [-R*ones(1,d);R*ones(1,d)];
[c4n,n4e] = distmeshnd(@s,@phi,h_min,box,fixed);

function val = s(x)
global r;
val = sqrt(sum(x.^2,2))-r;

function val = phi(x)
global alpha;
val = (1+sqrt(sum(x.^2,2))).^alpha;
