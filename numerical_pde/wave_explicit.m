function wave_explicit
T = 2; c = 1;
J = 20; K = 40;
Delta_x = 1/J;
Delta_t = T/K;
mu = c*Delta_t/Delta_x;
U = zeros(K+1,J+1);
for j = 0:J
    U(1,j+1) = u_0(j*Delta_x);
end
plot(Delta_x*(0:J),U(1,:)); axis([0,1,-1,1]); pause
U(2,1) = 0;
for j = 1:J-1
    U(2,j+1) = (1-mu^2)*U(1,j+1)+(mu^2/2)*(U(1,j)+U(1,j+2))...
            +Delta_t*v_0(j*Delta_x);
end
U(2,J+1) = 0;
for k = 1:K-1
    U(k+2,1) = 0;
    for j = 1:J-1
        U(k+2,j+1) = 2*(1-mu^2)*U(k+1,j+1)-U(k,j+1)...
            +mu^2*(U(k+1,j)+U(k+1,j+2));
    end
    U(k+2,J+1) = 0;
    plot(Delta_x*(0:J),U(k+2,:)); axis([0 1 -1 1]); pause(.1)
end
mesh(Delta_x*(0:J),Delta_t*(0:K),U);

function val = u_0(x)
% val = 2*max(min(x-1/4,3/4-x),0); 
val = sin(pi*x);

function val = v_0(x)
val = 0;
