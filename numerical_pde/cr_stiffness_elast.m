function [S_cr,b,vol_S] = cr_stiffness_elast(c4n,n4e,s4e,nS,d)
global lambda mu;
nE = size(n4e,1); ctr = 0; ctr_max = d^2*(d+1)^2*nE;
b = zeros(d*nS,1); vol_S = zeros(nS,1); zz = zeros(d,d);
I = zeros(ctr_max,1); J = zeros(ctr_max,1); X = zeros(ctr_max,1);
for j = 1:nE
    X_T = [ones(1,d+1);c4n(n4e(j,:),:)'];
    grads_T = X_T\[zeros(1,d);eye(d)];
    vol_T = det(X_T)/factorial(d);
    f_mp_T = f(sum(c4n(n4e(j,:),:),1)/(d+1));
    heights = 1./sqrt(sum(grads_T.^2,2));
    vol_S(s4e(j,:)) = factorial(d)*vol_T./heights;
    for m = 1:d+1
        for p = 1:d
            d_phi_mp = zz; d_phi_mp(p,:) = -d*grads_T(m,:);
            Ceps_phi_mp = lambda*sum(diag(d_phi_mp))*eye(d)...
                +2*mu*(d_phi_mp'+d_phi_mp)/2;
            b(d*(s4e(j,m)-1)+p) = b(d*(s4e(j,m)-1)+p)...
                +(1/(d+1))*vol_T*f_mp_T(p);
            for n = 1:d+1
                for q = 1:d
                   ctr = ctr+1; 
                   d_phi_nq = zz; d_phi_nq(q,:) = -d*grads_T(n,:);
                   I(ctr) = d*(s4e(j,m)-1)+p; 
                   J(ctr) = d*(s4e(j,n)-1)+q;
                   X(ctr) = vol_T*sum(sum(Ceps_phi_mp.*d_phi_nq));
                end
            end
        end
    end
end
S_cr = sparse(I,J,X,d*nS,d*nS);

function val = f(x); d = size(x,2); val = zeros(d,1); val(d) = -1;