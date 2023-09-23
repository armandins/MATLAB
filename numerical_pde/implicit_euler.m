function implicit_euler
T = 1;
J = 20; K = 10;
Delta_x = 1/J; Delta_t = T/K;
lambda = Delta_t/Delta_x^2;
U = zeros(K+1,J+1);
e = ones(J-1,1);
A = speye(J-1)+lambda*spdiags([-e,2*e,-e],[-1,0,1],J-1,J-1);
for j = 0:J
    U(1,j+1) = u_0(j*Delta_x);
end
plot(Delta_x*(0:J),U(1,:)); axis([0,1,0,1]); pause
for k = 0:K-1
    U(k+2,1) = 0;
    x = A\U(k+1,2:J)';
    U(k+2,2:J) = x';
    U(k+2,J+1) = 0;
    plot(Delta_x*(0:J),U(k+2,:)); axis([0,1,0,1]); pause(.05)
end
mesh(Delta_x*(0:J),Delta_t*(0:K),U);

function val = u_0(x)
val = sin(pi*x);