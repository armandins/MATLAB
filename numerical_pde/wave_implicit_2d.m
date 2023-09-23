function wave_implicit_2d
T = 100; c = 1;
J = 20; K = 20; L = (J-1)^2;
Delta_x = 1/J; Delta_t = T/K;
mu = c*Delta_t/Delta_x;
U = zeros(K+1,L);
e = ones(J-1,1); E = ones(L,1);
X = spdiags([-e,4*e,-e],[-1,0,1],J-1,J-1);
A = sparse(L,L);
for j = 1:(J-1):L
    A(j:j+J-2,j:j+J-2) = X;
end
A = A+spdiags([-E,-E],[-J+1,J-1],L,L);
B = speye(L)+(mu/4)*A;
C = speye(L)-(mu/4)*A;
for j = 1:J-1
    for m = 1:J-1
        U(1,j+(m-1)*(J-1)) = u_0([j,m]*Delta_x);
        V(j+(m-1)*(J-1)) = v_0([j,m]*Delta_x);
        F(j+(m-1)*(J-1)) = f(0,[j,m]*Delta_x);
    end
end
U(2,:) = B\(C*U(1,:)'+B*Delta_t*V'+(1/2)*Delta_t^2*F');
show(U(1,:),J,Delta_x); pause
for k = 1:K-1
    for j = 1:J-1
        for m = 1:J-1
            F((m-1)*(J-1)+j) = f(k*Delta_t,[j,m]*Delta_x);
        end
    end
    U(k+2,:) = B\(2*C*U(k+1,:)'-B*U(k,:)'+Delta_t^2*F');
    show(U(k+2,:),J,Delta_x); pause(.05);
end

function val = u_0(x)
val = 16*x(1)*(1-x(1))*x(2)*(1-x(2));

function val = v_0(x)
val = 0;

function val = f(t,x)
val = 0;

function show(U,J,Delta_x)
U_mat = zeros(J+1,J+1); U_mat(2:J,2:J) = reshape(U,J-1,J-1)';
mesh(Delta_x*(0:J),Delta_x*(0:J),U_mat);axis([0,1,0,1,-1,1,-1,1]); 
