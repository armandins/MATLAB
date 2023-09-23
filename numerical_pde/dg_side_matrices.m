function [s_sides,s_penal] = ...
    dg_side_matrices(n4e,e4s,n4s,vol_S,grads,normals_S)
global beta gamma;
nE = size(n4e,1); [nS,d] = size(n4s); 
ctr1 = 0; ctr2 = 0; ctr_max = 4*d^2*nS;
I = zeros(ctr_max,1); J = zeros(ctr_max,1); X = zeros(ctr_max,1);
K = zeros(ctr_max,1); L = zeros(ctr_max,1); Y = zeros(ctr_max,1);
shift = [2,1];
for j = 1:nS  
    h_S = vol_S(j)^(1/(d-1));
    for k = 1:2
        for ell = 1:2
            sigma_T_ell = 2*(e4s(j,ell)<e4s(j,shift(ell)))-1;
            if (e4s(j,k)&&e4s(j,ell))                
                [~,n] = setdiff(n4e(e4s(j,ell),:),n4s(j,:));
                for m = 1:d+1
                    ctr1 = ctr1+1;
                    I(ctr1) = (d+1)*(e4s(j,k)-1)+m;
                    J(ctr1) = (d+1)*(e4s(j,ell)-1)+n;
                    X(ctr1) = (-d*vol_S(j)/2)...
                        *grads((d+1)*(e4s(j,k)-1)+m,:)...
                        *normals_S(j,:)'*sigma_T_ell;
                end               
                [~,m] = setdiff(n4e(e4s(j,k),:),n4s(j,:));
                [~,n] = setdiff(n4e(e4s(j,ell),:),n4s(j,:));
                ctr2 = ctr2+1;
                K(ctr2) = (d+1)*(e4s(j,k)-1)+m; 
                L(ctr2) = (d+1)*(e4s(j,ell)-1)+n;
                Y(ctr2) = (-1)^(k-ell)*beta*h_S^(-gamma)*vol_S(j);
                [~,ind1] = intersect(n4e(e4s(j,k),:),n4s(j,:));
                [~,ind2] = intersect(n4e(e4s(j,ell),:),n4s(j,:));
                for m = reshape(ind1,1,d)
                   for n = reshape(ind2,1,d)
                     ctr2 = ctr2+1;
                     delta = (n4e(e4s(j,k),m)==n4e(e4s(j,ell),n));
                     K(ctr2) = (d+1)*(e4s(j,k)-1)+m; 
                     L(ctr2) = (d+1)*(e4s(j,ell)-1)+n;
                     Y(ctr2) = (-1)^(k-ell)*beta*h_S^(-gamma)...
                         *vol_S(j)*(d*delta-1)/(d+1);
                   end
                end              
            end
        end
    end
end
s_sides = sparse(I(1:ctr1),J(1:ctr1),X(1:ctr1),(d+1)*nE,(d+1)*nE);
s_penal = sparse(K(1:ctr2),L(1:ctr2),Y(1:ctr2),(d+1)*nE,(d+1)*nE);
