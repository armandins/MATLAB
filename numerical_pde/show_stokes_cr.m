function show_stokes_cr(c4n,n4e,s4e,u,p,mp_T)
[nC,d] = size(c4n);
u1_T = sum(u(d*(s4e-1)+1),2)/(d+1); 
u2_T = sum(u(d*(s4e-1)+2),2)/(d+1);
if d == 2
    trisurf(n4e,c4n(:,1),c4n(:,2),zeros(nC,1)',p); 
    hold on; view(0,90);
    quiver(mp_T(:,1),mp_T(:,2),u1_T,u2_T,'k');
    shading flat; hold off;
else
    tetramesh(n4e,c4n,p); hold on;
    u3_T = sum(u(d*s4e-1)+3,2)/(d+1);
    quiver3(mp_T(:,1),mp_T(:,2),mp_T(:,3),u1_T,u2_T,u3_T,'k');
    shading flat; hold off;
end
colorbar; 
