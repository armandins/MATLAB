function [edges,el2edges,Db2edges,Nb2edges] = edge_data(n4e,Db,Nb)
[nE,nV] = size(n4e); d = nV-1;
idx = [0,1,3,6]; Bdy = [Db;Nb]; nEdges = idx(d+1)*nE;
nDb = idx(d)*size(Db,1); nNb = idx(d)*size(Nb,1);
switch d
 case 1; edges = n4e;
 case 2; edges = [reshape(n4e(:,[1 2,1 3,2 3])',2,[])';Bdy];
 case 3; edges = reshape(n4e(:,[1 2,1 3,1 4,2 3,2 4,3 4])',2,[])';
     edges = [edges;reshape(Bdy(:,[1 2,1 3,2 3])',2,[])'];
end
[edges,~,edgeNumbers] = unique(sort(edges,2),'rows','first');
el2edges = reshape(edgeNumbers(1:nEdges),idx(d+1),[])';
Db2edges = reshape(edgeNumbers(nEdges+(1:nDb))',idx(d),[])';
Nb2edges = reshape(edgeNumbers(nEdges+nDb+(1:nNb))',idx(d),[])';
