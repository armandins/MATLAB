function theta_method_2d
T = 1;
theta = .5;
J = 20; K = 10; L = (J-1)^2;
Delta_x = 1/J; Delta_t = T/K;
lambda = Delta_t/Delta_x^2;
U = zeros(K+1,L);
e = ones(J-1,1); E = ones(L,1);
X = spdiags([-e,4*e,-e],[-1,0,1],J-1,J-1);
A = sparse(L,L);
for j = 1:(J-1):L
    A(j:j+J-2,j:j+J-2) = X;
end
A = A+spdiags([-E,-E],[-J+1,J-1],L,L);
B = speye(L)+lambda*theta*A;
C = speye(L)-lambda*(1-theta)*A;
for j = 1:J-1
    for m =1:J-1
        U(1,j+(m-1)*(J-1)) = u_0([j,m]*Delta_x);
    end
end
show(U(1,:),J,Delta_x); pause
for k = 0:K-1
    for j = 1:J-1
        for m = 1:J-1
            F(j+(m-1)*(J-1)) = f((k+theta)*Delta_t,[j,m]*Delta_x);
        end
    end
    U(k+2,:) = B\(C*U(k+1,:)'+Delta_t*F');
    show(U(k+2,:),J,Delta_x); pause(.05);
end

function val = u_0(x)
val = 16*x(1)*(1-x(1))*x(2)*(1-x(2));

function val = f(t,x)
val = 8;

function show(U,J,Delta_x)
U_mat = zeros(J+1,J+1); U_mat(2:J,2:J) = reshape(U,J-1,J-1)';
mesh(Delta_x*(0:J),Delta_x*(0:J),U_mat); axis([0,1,0,1,0,1,0,1]); 
