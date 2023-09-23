function show_p1(c4n,n4e,Db,Nb,u)
d = size(c4n,2);
if d == 1
    plot(c4n(n4e),u(n4e));
elseif d == 2
    trisurf(n4e,c4n(:,1),c4n(:,2),u);
elseif d == 3
    trisurf([Db;Nb],c4n(:,1),c4n(:,2),c4n(:,3),u);
end
