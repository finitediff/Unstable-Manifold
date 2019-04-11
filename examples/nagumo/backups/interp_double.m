function out = interp_double(sol,x,n)

temp = deval(sol,x);

mu_L = sqrt(3*n+1);
mu_R = -sqrt(3*n+1);

out = [exp(mu_L*x).*temp(1:2,:); exp(-mu_R*x).*temp(3:4,:)];



