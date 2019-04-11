function E = evans_co(s,p,lambda,e,m)





mu1_lam = sqrt(lambda+p.alpha);
mu2_lam = sqrt(lambda+1/p.gamma);

Q1 = [1;mu1_lam;0;0];
Q2 = [0;0;1;mu2_lam];


WL = [Q1, Q2];

OmegaL0 = orth(WL);
alphaL = OmegaL0'*WL;
muL = trace(OmegaL0'*e.LA(e.Li(1),lambda,s,p)*OmegaL0);
[omegal,gammal] = manifold_polar(e.Li,OmegaL0,lambda,e.LA,s,p,m,e.kl,muL);


Id = eye(4);
e1 = Id(:,1);
% e2 = Id(:,2);
e3 = Id(:,3);
% e4 = Id(:,4);


E = (det(alphaL)*gammal)*det([omegal e1,e3]);

