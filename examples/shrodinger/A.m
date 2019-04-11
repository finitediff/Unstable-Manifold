function out = A(x,lambda,s,p)

q = sqrt(2)*sech(x);

out = [ 0 0 1 0;
    0 0 0 1;
    1-3*q^2, lambda+2*p.nu, 0, 0;
    -lambda, p.mu-q^2, 0, 0];


