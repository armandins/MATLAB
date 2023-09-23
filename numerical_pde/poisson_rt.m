function poisson_rt(red)
[c4n,n4e,Db,Nb] = triang_cube(2);
for j = 1:red
    [c4n,n4e,Db,Nb] = red_refine(c4n,n4e,Db,Nb);
end
[s4e,sign_s4e,n4s,s4Db,s4Nb] = sides(n4e,Db,Nb);
nS = size(n4s,1); nE = size(n4e,1); 
nNb = size(Nb,1); nDb = size(Db,1); 
m_loc = [2,1,1;1,2,1;1,1,2]/12;
A = sparse(nS,nS); B = sparse(nS,nE); D = sparse(nS,nNb);
b = zeros(nS+nE+nNb,1); 
shift1 = [2,3,1]; shift2 = [3,1,2];
for j = 1:nE
    area_T = det([1,1,1;c4n(n4e(j,:),:)'])/2;
    mp_T = sum(c4n(n4e(j,:),:))/3;
    for m = 1:3
        length_m = norm(c4n(n4e(j,shift1(m)),:)...
           -c4n(n4e(j,shift2(m)),:));
        for n = 1:3
            length_n = norm(c4n(n4e(j,shift1(n)),:)...
                -c4n(n4e(j,shift2(n)),:));
            for o = 1:3
                for p = 1:3
                    A(s4e(j,m),s4e(j,n)) = A(s4e(j,m),s4e(j,n))...
                        +length_m*length_n/(4*area_T)...
                        *sign_s4e(j,m)*sign_s4e(j,n)*m_loc(o,p)...
                        *(c4n(n4e(j,m),:)-c4n(n4e(j,o),:))...
                        *(c4n(n4e(j,n),:)-c4n(n4e(j,p),:))';
                end
            end
        end
        B(s4e(j,m),j) = -length_m*sign_s4e(j,m);
        b(nS+j) = -area_T*g(mp_T);
    end
end
for j = 1:nDb
    length_S = norm(c4n(Db(j,2),:)-c4n(Db(j,1),:))';
    mp_S = (c4n(Db(j,1),:)+c4n(Db(j,2),:))/2;    
    b(s4Db(j)) = p_D(mp_S)*length_S;
end
for j = 1:nNb
    mp_S = (c4n(Nb(j,1),:)+c4n(Nb(j,2),:))/2;
    D(s4Nb(j),j) = 1; b(nS+nE+j) = sigma(mp_S);
end
O1 = sparse(nE,nE+nNb); O2 = sparse(nNb,nE+nNb);
G = [A,B,D;B',O1;D',O2];
x = G\b; u = x(1:nS); p = x(nS+(1:nE));
show_rt(c4n,n4e,u,p,s4e,sign_s4e);

function val = sigma(x); val = ones(size(x,1),1);
function val = p_D(x); val = sin(2*pi*x(:,1));
function val = g(x); val = ones(size(x,1),1);
