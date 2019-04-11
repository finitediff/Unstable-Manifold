function out = bc(ya,yb,n,lam,p)

% We get eigenvalues for A+ = A-,
% lambda = -sqrt(μ - sqrt(μ^2 - 2*μ + 4*λ^2*n^2 + 8*λ*ν*n + 1) + 1)/sqrt(2)
% lambda = -sqrt(μ + sqrt(μ^2 - 2*μ + 4*λ^2*n^2 + 8*λ*ν*n + 1) + 1)/sqrt(2)
% lambda = sqrt(μ - sqrt(μ^2 - 2*μ + 4*λ^2*n^2 + 8*λ*ν*n + 1) + 1)/sqrt(2)
% lambda = sqrt(μ + sqrt(μ^2 - 2*μ + 4*λ^2*n^2 + 8*λ*ν*n + 1) + 1)/sqrt(2)

mu1 = sqrt(p.mu - sqrt(p.mu^2 - 2*p.mu + 4*lam^2*n^2 + 8*lam*p.nu*n + 1) + 1)/sqrt(2);
mu2 = sqrt(p.mu + sqrt(p.mu^2 - 2*p.mu + 4*lam^2*n^2 + 8*lam*p.nu*n + 1) + 1)/sqrt(2);

PL = [mu1,-1,0,0;
      0, 0, mu2, -1];
  
PR = [-mu1,-1,0,0;
      0, 0, -mu2,-1];

% out = [PL*ya;yb(2);yb(4)];

% out = [PL*ya;yb(1);yb(3)];


out = [PL*ya;PR*yb];











