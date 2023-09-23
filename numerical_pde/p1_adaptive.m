function p1_adaptive(red) % Nb = []; d = 2; 
theta = .5; eps_stop = 5e-2; error_bound = 1;
c4n = [-1,-1;0,-1;-1,0;0,0;1,0;-1,1;0,1;1,1];
n4e = [1,2,4;4,3,1;3,4,7;7,6,3;4,5,8;8,7,4];
Db = [1,2;2,4;4,5;5,8;8,7;7,6;6,3;3,1]; 
for j = 1:red
    marked = ones(size(n4e,1),1);
    [c4n,n4e,Db] = rgb_refine(c4n,n4e,Db,marked);
end
while error_bound > eps_stop
    %%% solve
    fNodes = setdiff(1:size(c4n,1),unique(Db));
    u = zeros(size(c4n,1),1);
    [s,m] = fe_matrices(c4n,n4e);
    b = m*f(c4n);
    u(fNodes) = s(fNodes,fNodes)\b(fNodes);
    show_p1(c4n,n4e,Db,[],u); view(0,90); pause(.05)
    %%% estimate
    eta = comp_estimators(c4n,n4e,Db,u);
    error_bound = sqrt(sum(eta.^2))
    %%% mark
    marked = (eta>theta*max(eta));
    %%% refine
    if error_bound > eps_stop 
       [c4n,n4e,Db] = rgb_refine(c4n,n4e,Db,marked);
    end
end

function eta = comp_estimators(c4n,n4e,Db,u)
[s4e,~,n4s,s4Db] = sides(n4e,Db,[]);
[~,d] = size(c4n); nE = size(n4e,1); nS = size(n4s,1);
eta_S = zeros(nS,1); eta_T_sq = zeros(nE,1);
for j = 1:nE
    X_T = [ones(1,d+1);c4n(n4e(j,:),:)'];
    grads_T = X_T\[zeros(1,d);eye(d)];
    vol_T = det(X_T)/factorial(d);
    h_T = vol_T^(1/d);
    mp_T = sum(c4n(n4e(j,:),:),1)/(d+1);
    eta_T_sq(j) = h_T^2*vol_T*(f(mp_T));
    nabla_u_T = grads_T'*u(n4e(j,:));
    normal_times_area = -grads_T*vol_T*d;    
    eta_S(s4e(j,:)) = h_T^((2-d)/2)*eta_S(s4e(j,:))...
        +normal_times_area*nabla_u_T;
end
eta_S(s4Db) = 0;
eta_S_T_sq = sum(eta_S(s4e).^2,2);
eta = (eta_T_sq+eta_S_T_sq).^(1/2);
 
function val = f(x); val = ones(size(x,1),1);