function [c4n,n4e,Db] = rgb_refine(c4n,n4e,Db,marked)
% edges n4e(:,[1,3]) have to be longest edges 
[edges,el2edges,Db2edges] = edge_data_2d(n4e,Db);
nC = size(c4n,1); nEdges = size(edges,1); tmp = 1;
markedEdges = zeros(nEdges,1);
markedEdges(reshape(el2edges(marked==1,[1 2 3]),[],1)) = 1;
while tmp > 0
    tmp = nnz(markedEdges);
    el2markedEdges = markedEdges(el2edges);
    el2markedEdges(el2markedEdges(:,1)+el2markedEdges(:,3)>0,2)=1;
    markedEdges(el2edges(el2markedEdges==1))= 1;
    tmp = nnz(markedEdges)-tmp;
end
newNodes = zeros(nEdges,1);
newNodes(markedEdges==1) = (1:nnz(markedEdges))'+nC;
newInd = newNodes(el2edges);
red   = newInd(:,1) >  0 & newInd(:,2) >  0 & newInd(:,3) >  0;
blue1 = newInd(:,1) >  0 & newInd(:,2) >  0 & newInd(:,3) == 0;
blue3 = newInd(:,1) == 0 & newInd(:,2) >  0 & newInd(:,3) >  0;
green = newInd(:,1) == 0 & newInd(:,2) >  0 & newInd(:,3) == 0;
remain= newInd(:,1) == 0 & newInd(:,2) == 0 & newInd(:,3) == 0;
n4e_red = [n4e(red,1),newInd(red,[3 2]),...
           newInd(red,[2 1]),n4e(red,3),...
           newInd(red,3),n4e(red,2),newInd(red,1),...
           newInd(red,:)];   
n4e_red = reshape(n4e_red',3,[])';
n4e_blue1 = [n4e(blue1,2),newInd(blue1,2),n4e(blue1,1) ...
             n4e(blue1,2),newInd(blue1,[1 2]),...
             newInd(blue1,[2 1]),n4e(blue1,3)];
n4e_blue1 = reshape(n4e_blue1',3,[])';
n4e_blue3 = [n4e(blue3,1),newInd(blue3,3),newInd(blue3,2),...
             newInd(blue3,2),newInd(blue3,3),n4e(blue3,2),...
             n4e(blue3,3),newInd(blue3,2),n4e(blue3,2)];
n4e_blue3 = reshape(n4e_blue3',3,[])';
n4e_green = [n4e(green,2),newInd(green,2),n4e(green,1),...
             n4e(green,3),newInd(green,2),n4e(green,2)];
n4e_green = reshape(n4e_green',3,[])';
n4e = [n4e(remain,:);n4e_red;n4e_blue1;n4e_blue3;n4e_green];
newCoord =.5*(c4n(edges(markedEdges==1,1),:)...
    +c4n(edges(markedEdges==1,2),:));
c4n = [c4n;newCoord];
newDb = newNodes(Db2edges); ref = newDb>0; Db_old = Db(~ref,:);
Db_new = [Db(ref,1),newDb(ref),newDb(ref),Db(ref,2)];
Db_new = reshape(Db_new',2,[])'; Db = [Db_old;Db_new];

function [edges,el2edges,Db2edges] = edge_data_2d(n4e,Db)
nE = size(n4e,1); nEdges = 3*nE; nDb = size(Db,1);
edges = [reshape(n4e(:,[2 3,1 3,1 2])',2,[])';Db];
[edges,~,edgeNumbers] = unique(sort(edges,2),'rows','first');
el2edges = reshape(edgeNumbers(1:nEdges),3,[])';
Db2edges = reshape(edgeNumbers(nEdges+(1:nDb))',1,[])';
