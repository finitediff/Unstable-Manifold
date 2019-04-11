function out = bc_eigfun(ya,yb,lambda,p)





% mu1 = sqrt(lambda+p.alpha);
% mu2 = sqrt(lambda+1/p.gamma);
% 
% z1 = [-mu1;1;0;0];
% z2 = [0;0;-mu2;1];
% 
% z3 = [mu1;1;0;0];
% z4 = [0;0;mu2;1];
% 
% AM = s.Flinear(s.UL,p);
AM = A(-p.L);
PL = orth(projection1(AM,-1,0).').';
% AP = s.Flinear(s.UR,p);
AP = A(p.L);
PR = orth(projection1(AP,1,0).').';


out = [PL.'*ya(1:7);PR.'*yb(8:14);norm(yb(1:7))^2-1;yb(1:7)-ya(8:14)];



