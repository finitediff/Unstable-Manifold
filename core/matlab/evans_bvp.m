function E = evans_bvp(s,p,N,P,lambda,sig1,sig2,L)





mu1_lam = sqrt(lambda+p.alpha);
mu2_lam = sqrt(lambda+1/p.gamma);

Q1 = [1;mu1_lam;0;0];
Q2 = [0;0;1;mu2_lam];

Proj = null(Q2.').';

degree = 30;
m.sys_dim = 4;
m.stats = 'off';
mesh = [L,0.5*L,0];
PR = [];
PL = Proj;
mu = mu2_lam;
zj = Q2';Amat = @(x)(Av(x,lambda,s,p));
[vL2,vR2] = solve_bvp_cheb(m,Amat,mu,PR,PL,zj,mesh,degree);

% scl2 = sign(Q2(4))*sign(vL2(4))*norm(Q2)/norm(vL2)

ind = find(abs(Q2)== max(abs(Q2)));
scl2 = Q2(ind(1))/vL2(ind(1));

W2 = scl2*vR2;


PR = vR2';
PL = null([Q1,Q2].').';
mu = mu1_lam;
zj = Q1';
Amat = @(x)(Av(x,lambda,s,p));
[vL1,vR1] = solve_bvp_cheb(m,Amat,mu,PR,PL,zj,mesh,degree);

ind = find(abs(Q1)== max(abs(Q1)));
scl1 = Q1(ind(1))/vL1(ind(1));

W1 = scl1*vR1;

Id = eye(4);
e1 = Id(:,1);
e2 = Id(:,2);
e3 = Id(:,3);
e4 = Id(:,4);

E = det([W1,W2,e1,e3]);








