function [c4n,n4e,Db,Nb,P0,P1] = red_refine(c4n,n4e,Db,Nb)
nE = size(n4e,1); [nC,d] = size(c4n);
[edges,el2edges,Db2edges,Nb2edges] = edge_data(n4e,Db,Nb);
nS = size(edges,1);
newNodes = nC+(1:nS)';
newIndices = reshape(newNodes(el2edges),size(el2edges));
idxElements = refinement_rule(d);
n4e = [n4e,newIndices];
n4e = n4e(:,idxElements);
n4e = reshape(n4e',d+1,[])';
newCoord = .5*(c4n(edges(:,1),:)+c4n(edges(:,2),:));
c4n = [c4n;newCoord];
idx_Bdy = refinement_rule(d-1);
if size(Db,1) > 0
    newDb = reshape(newNodes(Db2edges),size(Db2edges));
    Db = [Db,newDb]; Db = Db(:,idx_Bdy); Db = reshape(Db',d,[])';
end
if size(Nb,1) > 0
    newNb = reshape(newNodes(Nb2edges),size(Nb2edges));
    Nb = [Nb,newNb]; Nb = Nb(:,idx_Bdy); Nb = reshape(Nb',d,[])';
end
P0 = sparse(1:size(n4e,1),reshape( repmat(1:nE,2^d,1),1,[]),1);
P1 = sparse([1:nC,newNodes',newNodes'],...
    [1:nC,edges(:,1)',edges(:,2)'],...
    [ones(1,nC),.5*ones(1,2*nS)],nC+nS,nC);

function idx = refinement_rule(d)
switch d
 case 0; idx = 1;
 case 1; idx = [1 3,3 2];
 case 2; idx = [1 4 5,5 6 3,4 2 6,6 5 4];
 case 3; idx = [1 5 6 7,5 2 8 9,6 8 3 10,7 9 10 4,...
            5 6 7 9,9 6 8 5,6 7 9 10,10 8 9 6];
end

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

