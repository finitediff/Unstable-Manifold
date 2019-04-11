function out = bc_eigfun(ya,yb,lambda,p,s)

% lambda

% A(s.L,lambda,s,p)


% [PL,QL] = projection2(A(s.L,lambda,s,p),1,0);
% [PR,QR] = projection2(A(s.L,lambda,s,p),-1,0);


Q = 0.5*sqrt((p.mu+1)^2-4*(lambda^2+2*lambda*p.nu+p.mu));
C = 0.5*(1+p.mu);
mu1 = (sqrt(Q+C));
mu2 = (sqrt(C-Q));

PL = [mu1,-1,0,0;
      0, 0, mu2, -1];
  
PR = [-mu1,-1,0,0;
      0, 0, -mu2,-1];


% PL = [mu1L, (mu1L^2-mu1L)/(conj(lambda)+2*p.nu), 1, (1-mu1L^2)/conj(lambda);
%     mu2L, (mu2L^2-mu2L)/(conj(lambda)+2*p.nu), 1, (1-mu2L^2)/conj(lambda)];
% 
% PR = [-mu1L, (mu1L^2+mu1L)/(conj(lambda)+2*p.nu), 1, (1-mu1L^2)/conj(lambda);
%     -mu2L, (mu2L^2+mu2L)/(conj(lambda)+2*p.nu), 1, (1-mu2L^2)/conj(lambda)];

  
out = [PL*ya(1:4);PR*yb(5:8);norm(yb(1:4))^2-1;yb(1:4)-ya(5:8)];

% out = [QL.'*ya(1:4);QR.'*yb(5:8);norm(yb(1:4))^2-1;yb(1:4)-ya(5:8)];



