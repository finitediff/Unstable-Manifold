function out = bc_fun(ya,yb,lambda,p)





mu1 = -sqrt(lambda+p.alpha);
mu2 = -sqrt(lambda+1/p.gamma);

z1 = [-mu1;1;0;0];
z2 = [0;0;-mu2;1];

e2 = [0;1;0;0];
e4 = [0;0;0;1];

out = [[z1,z2].'*ya;[e2,e4].'*yb;norm(yb)^2-1];



