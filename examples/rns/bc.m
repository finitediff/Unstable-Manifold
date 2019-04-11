%bc_L_fun = @(U_n,U_o,K,H,p)[U_n(1,1)-bc_L(1);U_n(2,1)-bc_L(2); ...
%    U_n(3,1)-bc_L(3);U_n(4,1)-bc_L(4)];
%bc_R_fun = @(U_n,U_o,K,H,p)[U_n(1,end)-bc_R(1);U_n(2,end)-bc_R(2); ...
%    U_n(3,end)-bc_R(3);U_n(4,end)-U_n(4,end-1)];
%bc_L_jac_fun = @(U_n,U_o,K,H,p)[1,0,0,0; 0,1,0,0; 0,0,1,0; 0,0,0,1];
%bc_R_jac_fun = @(U_n,U_o,K,H,p)[0 0 0 0 1,0,0,0; 0 0 0 0 0,1,0,0;...
% 0 0 0 0 0,0,1,0; 0 0 0 -1 0,0,0,1];



function out = bc(ya,yb,n,lam,p)

bc_L = [ya - p.L]
bc_R = [yb - p.R]

out = [bc_L, bc_R];











