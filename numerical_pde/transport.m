function transport
a = 1; T = 1;
J = 20; K = 20;
Delta_x = 1/J; Delta_t = T/K;
U = zeros(K+1,J+1);
for j = 0:J
    U(1,j+1) = u_0(j*Delta_x);
end
plot(Delta_x*(0:J),U(1,:)); axis([0,1,0,1]); pause
for k = 0:K-1
    U(k+1,1) = 0;
    for j = 1:J
        U(k+2,j+1) = U(k+1,j+1)...
            -a*(Delta_t/Delta_x)*(U(k+1,j+1)-U(k+1,j));
    end
    plot(Delta_x*(0:J),U(k+2,:)); axis([0,1,0,1]); pause(.1)
end
clf; mesh(Delta_x*(0:J),Delta_t*(0:K),U);

function val = u_0(x)
if x <= .5
    val = .5*(1+sin(4*pi*x-pi/2));
else
    val = 0;
end
