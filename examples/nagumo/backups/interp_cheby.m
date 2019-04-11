function out = interp_cheby(cf,x,n,N,a,b)


x0 = (x - 0.5*(a+b))/(0.5*(a-b));

theta = acos(x0);

T = cos(theta.'*(0:1:N-1));
temp = (T*cf).';


mu_L = sqrt(3*n+1);
mu_R = -sqrt(3*n+1);

out = [exp(mu_L*x).*temp(1:2,:); exp(-mu_R*x).*temp(3:4,:)];



