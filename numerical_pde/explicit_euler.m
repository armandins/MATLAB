function explicit_euler
T = 1;
J = 20; K = 400;
Delta_x = 1/J; Delta_t = T/K;
lambda = Delta_t/Delta_x^2;
U = zeros(K+1,J+1);
for j = 0:J
    U(1,j+1) = u_0(j*Delta_x);
end
plot(Delta_x*(0:J),U(1,:)); axis([0,1,0,1]); pause
for k = 0:K-1
    U(k+2,1) = 0;
    for j = 1:J-1
        U(k+2,j+1) = (1-2*lambda)*U(k+1,j+1)...
            +lambda*(U(k+1,j+2)+U(k+1,j));
    end
    U(k+1,J+1) = 0;
    plot(Delta_x*(0:J),U(k+2,:)); axis([0,1,0,1]); pause(.1)
end
clf; mesh(Delta_x*(0:J),Delta_t*(0:K),U);

function val = u_0(x)
val = sin(pi*x);
