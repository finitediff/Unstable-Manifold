function out = Av(x,lambda,s,p)

% Evans function matrix
%
% W'(x) = A(x,\lambda)W(x)


temp = deval(s.sol,-abs(x));
n = size(temp,2);
u = temp(1,:);
v = temp(3,:);

out = zeros(4,4,n);

one = ones(1,n);
out(1,2,:) = one;
out(2,1,:) = lambda+v.^2+p.alpha;
out(2,3,:) = 2*u.*v;
out(3,4,:) = one;
out(4,1,:) = -v.^2/p.gamma;
out(4,3,:) = lambda+(1-2*u.*v)/p.gamma;

% out = [0 1 0 0; 
%     lambda+v^2+p.alpha 0 2*u*v 0;
%     0 0 0 1;
%     -v^2/p.gamma 0 lambda+(1-2*u*v)/p.gamma 0];



