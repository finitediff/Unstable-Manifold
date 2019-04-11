function out = profile_guess(x,p)

% use explicit solution as an initial guess
xi = x/sqrt(p.gamma);
if 1-9*p.gamma/2 > 0
    Q = sqrt(1-9*p.gamma/2);
else
    Q = 0.05;
end


out = [1-3*p.gamma/(1+Q*cosh(xi)); 
    3*Q*sqrt(p.gamma)*sinh(xi)/(1+Q*cosh(xi))^2;
    3/(1+Q*cosh(xi));
    -3*Q*sinh(xi)*(1/sqrt(p.gamma))/(1+Q*cosh(xi))^2];

% hold on
% plot(x,out,'.k')









