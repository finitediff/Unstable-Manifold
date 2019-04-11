function E = evans_parametrized(s,p,N,P,lambda,sig1,sig2,L)




[V1,V2,err1,err2] = evans_parametrization(p,N,P,lambda);
 
% fprintf('\n\nMax residual error for parametrization of Evans function: %4.4g, %4.4g\n\n',...
%     err1,err2);


Q1 = eval_P(V1,N,N,sig1,sig2);

Q2 = eval_P(V2,N,N,sig1,sig2);

% toc(t1)

mu1_lam = sqrt(lambda+p.alpha);
mu2_lam = sqrt(lambda+1/p.gamma);



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








