function R = high_frequency_bound(p,s)


x = linspace(s.L,0,1000);
u = zeros(size(x));
v = u;

for j = 1:length(x)
    temp = deval(s.sol,x(j));
    u(j) = temp(1);
    v(j) = temp(3);
end

u_inf = 0;
v_inf = 0;
for j = 1:length(x)
   u_inf = max(u_inf,abs(u(j)));
   v_inf = max(v_inf,abs(v(j)));
end

Re_lambda = max(0, max(v_inf^2+2*u_inf*v_inf-p.alpha, ...
    (1/abs(p.gamma))*v_inf^2+(2/abs(p.gamma))*u_inf*v_inf-1/p.gamma));

Im_lambda = sqrt( 2*v_inf^3*u_inf/abs(p.gamma));

R = sqrt(Re_lambda^2+Im_lambda^2);









