function [s4e,sign_s4e,n4s,s4Db,s4Nb,e4s] = sides(n4e,Db,Nb)
nE = size(n4e,1); d = size(n4e,2)-1; 
nDb = size(Db,1); nNb = size(Nb,1); 
if d == 2
    Tsides = [n4e(:,[2,3]);n4e(:,[3,1]);n4e(:,[1,2])];
else
    Tsides = [n4e(:,[2,4,3]);n4e(:,[1,3,4]);...
        n4e(:,[1,4,2]);n4e(:,[1,2,3])];
end
[n4s,i2,j] = unique(sort(Tsides,2),'rows');
s4e = reshape(j,nE,d+1); nS = size(n4s,1);
sign_s4e = ones((d+1)*nE,1); sign_s4e(i2) = -1; 
sign_s4e = reshape(sign_s4e,nE,d+1);
[~,~,j2] = unique(sort([n4s;Db;Nb],2),'rows'); 
s4Db = j2(nS+(1:nDb)); s4Nb = j2(nS+nDb+(1:nNb));
e4s = zeros(size(n4s,1),2);
e4s(:,1) = mod(i2-1,nE)+1;
i_inner = setdiff(1:(d+1)*nE,i2);
e4s(j(i_inner),2) = mod(i_inner-1,nE)+1;