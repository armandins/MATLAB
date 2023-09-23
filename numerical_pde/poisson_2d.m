function poisson_2d
J = 20;
Delta_x = 1/J;
L = (J-1)^2;
e = ones(J-1,1); E = ones(L,1);
X = spdiags([-e,4*e,-e],[-1,0,1],J-1,J-1);
A = sparse(L,L);
for j = 1:(J-1):L
    A(j:j+J-2,j:j+J-2) = X;
end
A = A+spdiags([-E,-E],[-J+1,J-1],L,L);
for j = 1:J-1
    for m = 1:J-1
        b(j+(m-1)*(J-1),1) = Delta_x^2*f([j,m]*Delta_x);
    end
end
U = A\b;
show(U,J,Delta_x);

function val = f(x)
val = 1;

function show(U,J,Delta_x)
U_mat = zeros(J+1,J+1); U_mat(2:J,2:J) = reshape(U,J-1,J-1)';
mesh(Delta_x*(0:J),Delta_x*(0:J),U_mat); axis([0,1,0,1,0,.1]); 