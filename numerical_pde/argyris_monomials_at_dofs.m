function N = argyris_monomials_at_dofs(loc_c4n)
J = [0,-1;1,0];
x = loc_c4n(:,1); y = loc_c4n(:,2);
e = ones(3,1); o = zeros(3,1);
p = [e,x,y,x.^2,x.*y,y.^2,x.^3,x.^2.*y,x.*y.^2,...
    y.^3,x.^4,x.^3.*y,x.^2.*y.^2,x.*y.^3,y.^4,...
    x.^5,x.^4.*y,x.^3.*y.^2,x.^2.*y.^3,x.*y.^4,y.^5];
p_x = [o,e,o,2*x,y,o,3*x.^2,2*x.*y,y.^2,o,...
    4*x.^3,3*x.^2.*y,2*x.*y.^2,y.^3,o,...
    5*x.^4,4*x.^3.*y,3*x.^2.*y.^2,2*x.*y.^3,y.^4,o];
p_y = [o,o,e,o,x,2*y,o,x.^2,2*x.*y,3*y.^2,...
    o,x.^3,2*x.^2.*y,3*x.*y.^2,4*y.^3,...
    o,x.^4,2*x.^3.*y,3*x.^2.*y.^2,4*x.*y.^3,5*y.^4];
p_xx = [o,o,o,2*e,o,o,6*x,2*y,o,o,...
    12*x.^2,6*x.*y,2*y.^2,o,o,...
    20*x.^3,12*x.^2.*y,6*x.*y.^2,2*y.^3,o,o];
p_yy = [o,o,o,o,o,2*e,o,o,2*x,6*y,...
    o,o,2*x.^2,6*x.*y,12*y.^2,...
    o,o,2*x.^3,6*x.^2.*y,12*x.*y.^2,20*y.^3];
p_xy = [o,o,o,o,e,o,o,2*x,2*y,o,o,3*x.^2,4*x.*y,3*y.^2,o,...
    o,4*x.^3,6*x.^2.*y,6*x.*y.^2,4*y.^3,o];
N = [p;p_x;p_y;p_xx;p_yy;p_xy];
shift1 = [2,3,1]; shift2 = [3,1,2];
for k = 1:3
    z_a = loc_c4n(shift2(k),:); z_b = loc_c4n(shift1(k),:);      
    m_S = (z_a+z_b)/2; n_S = J*(z_b-z_a)'/norm(z_b-z_a);
    x = m_S(1); y = m_S(2); e = 1; o = 0; 
    p_x = [o,e,o,2*x,y,o,3*x.^2,2*x.*y,y.^2,o,...
        4*x.^3,3*x.^2.*y,2*x.*y.^2,y.^3,o,...
        5*x.^4,4*x.^3.*y,3*x.^2.*y.^2,2*x.*y.^3,y.^4,o];
    p_y = [o,o,e,o,x,2*y,o,x.^2,2*x.*y,3*y.^2,...
        o,x.^3,2*x.^2.*y,3*x.*y.^2,4*y.^3,...
        o,x.^4,2*x.^3.*y,3*x.^2.*y.^2,4*x.*y.^3,5*y.^4];
    N(18+k,:) = n_S(1)*p_x+n_S(2)*p_y;
end
