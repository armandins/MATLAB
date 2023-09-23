function show_chorin(c4n,n4e,u,p)
d = size(c4n,2);
if d == 2
    trisurf(n4e,c4n(:,1),c4n(:,2),p-max(p)-1); hold on;
    quiver(c4n(:,1),c4n(:,2),u(1:2:end),u(2:2:end),'k.'); 
    hold off; view(0,90); shading flat; drawnow;
elseif d == 3
    p_T = sum(p(n4e),2)/(d+1);
    idx = find(c4n(n4e(:,1),2)>0);
    tetramesh(n4e(idx,:),c4n,p_T(idx)); hold on;
    quiver3(c4n(:,1),c4n(:,2),c4n(:,3),...
        u(1:3:end),u(2:3:end),u(3:3:end),'k'); 
    view(-26,14); hold off; drawnow;
end
