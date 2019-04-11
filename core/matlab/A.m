function out = A(x,lambda,s,p)

% Evans function matrix
%
% W'(x) = A(x,\lambda)W(x)


temp = deval(s.sol,-abs(x));
u = temp(1);
v = temp(3);

out = [0 1 0 0; 
    lambda+v^2+p.alpha 0 2*u*v 0;
    0 0 0 1;
    -v^2/p.gamma 0 lambda+(1-2*u*v)/p.gamma 0];



