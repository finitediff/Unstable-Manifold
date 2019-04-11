function out = profile_jacobian(y,p)

v1 = y(1);
v2 = y(3);

out = [ 0 1 0 0;
    p.L^2*(v2^2+p.alpha) 0 2*p.L^2*v1*v2 0;
    0 0 0 1;
    -p.L^2*v2^2/p.gamma 0 p.L^2*(1-2*v1*v2)/p.gamma 0];








