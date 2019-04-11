function out = bc_parametrization(ya,yb,var,P,N)


sig1 = var(1);
sig2 = var(2);

z = eval_P(P,N,N,sig1,sig2);

out = [ya-z;yb(2);yb(4)];