function dirichlet_neumann(red)
c4n{1} = [0 0;1 0;0 1;1 1]; 
c4n{2} = [1 0;2 0;1 1;2 1]; 
n4e{1} = [1 2 4;1 4 3]; 
n4e{2} = [1 2 4;1 4 3]; 
Db{1} = [1 2;4 3;3 1]; 
Db{2} = [1 2;2 4;4 3];
for j = 1:red
    [c4n{1},n4e{1},Db{1}] = red_refine(c4n{1},n4e{1},Db{1},[]);
    [c4n{2},n4e{2},Db{2}] = red_refine(c4n{2},n4e{2},Db{2},[]);
end
for j = 1:2
    nC{j} = size(c4n{j},1);
    dNodes{j} = unique(Db{j});
    [s{j},m{j}] = fe_matrices(c4n{j},n4e{j});
    b{j} = m{j}*f(c4n{j});
end
[~,gamma{1},gamma{2}] = intersect(c4n{1},c4n{2},'rows');
dNodes{1} = union(dNodes{1},gamma{1});
for j = 1:2
    fNodes{j} = setdiff(1:nC{j},dNodes{j});
end
lambda = zeros(length(gamma{1}),1);
theta = 1/4; 
eps_stop = 1e-3; diff = 1; 
while diff > eps_stop
    %%% Initialize 
    u{1} = zeros(nC{1},1); u{2} = zeros(nC{2},1);    
    %%% Step (1) 
    u{1}(gamma{1}) = lambda;
    b1 = b{1}-s{1}*u{1};
    u{1}(fNodes{1}) = s{1}(fNodes{1},fNodes{1})\b1(fNodes{1});    
    %%% Step (2) 
    b2 = b{2};
    normal_trans = b{1}-s{1}*u{1};
    b2(gamma{2}) = b2(gamma{2})+normal_trans(gamma{1});
    u{2}(fNodes{2}) = s{2}(fNodes{2},fNodes{2})\b2(fNodes{2});    
    %%% Step (3) 
    lambda = theta*u{2}(gamma{2})+(1-theta)*lambda;
    diff = max(abs((u{1}(gamma{1})-u{2}(gamma{2}))))
    %%% Visualize
    show_p1(c4n{1},n4e{1},Db{1},[],u{1}); hold on; 
    show_p1(c4n{2},n4e{2},Db{2},[],u{2}); hold off; 
    pause(1);
end

function val = f(x); val = ones(size(x,1),1);

