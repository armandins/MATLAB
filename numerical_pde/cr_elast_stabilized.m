function cr_elast_stabilized(d,red)
global lambda mu;
lambda = 1000; mu = 1; beta = 2*mu;
[c4n,n4e,Db,Nb] = triang_cube(d);
for j = 1:red
    [c4n,n4e,Db,Nb] = red_refine(c4n,n4e,Db,Nb);
end
[s4e,~,n4s,s4Db,~,e4s] = sides(n4e,Db,Nb);
nS = size(n4s,1); 
fSides = setdiff(1:nS,s4Db);
FSides = d*(repmat(fSides,d,1)-1)+(1:d)'*ones(1,size(fSides,2));
FSides = FSides(:);
u = zeros(d*nS,1); 
[S_cr,b,vol_S] = cr_stiffness_elast(c4n,n4e,s4e,nS,d);
ctr = 0; ctr_max = 4*d^3*nS;
I = zeros(ctr_max,1); J = zeros(ctr_max,1); X = zeros(ctr_max,1);
for j = 1:nS  
    h_S = vol_S(j)^(1/(d-1));
    if (e4s(j,1) && e4s(j,2))
        for k = 1:2
            for ell = 1:2
                [~,ind1] = intersect(n4e(e4s(j,k),:),n4s(j,:));
                [~,ind2] = intersect(n4e(e4s(j,ell),:),n4s(j,:));
                for m = reshape(ind1,1,d)
                    for n = reshape(ind2,1,d)
                        for p = 1:d
                            ctr = ctr+1;
                            delta = (n4e(e4s(j,k),m)...
                                ==n4e(e4s(j,ell),n));
                            I(ctr) = d*(s4e(e4s(j,k),m)-1)+p;
                            J(ctr) = d*(s4e(e4s(j,ell),n)-1)+p;
                            X(ctr) = (-1)^(k-ell)*beta*h_S^(-1)...
                                *vol_S(j)*(d*delta-1)/(d+1);
                        end
                    end
                end
            end
        end
    end
end
Z_stab = sparse(I(1:ctr),J(1:ctr),X(1:ctr),d*nS,d*nS);
S_full = S_cr+Z_stab;
u(FSides) = S_full(FSides,FSides)\b(FSides);
show_cr_def(c4n,n4e,s4e,u)
