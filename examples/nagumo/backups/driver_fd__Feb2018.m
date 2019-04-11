clc; close all; beep off;


dom = linspace(-10,10,101);

bc_L = sqrt(2)*sech(dom(1));
bc_R = sqrt(2)*sech(dom(end));

bc_L_fun = @(U_n,U_o,K,H,p)(U_n(1,1)-bc_L);
bc_R_fun = @(U_n,U_o,K,H,p)(U_n(1,end)-bc_R);
bc_L_jac_fun = @(U_n,U_o,K,H,p)1 ;
bc_R_jac_fun = @(U_n,U_o,K,H,p)1;

% bc_L_fun = @(U_n,U_o,K,H,p)(U_n(1,1)-U_n(1,2));
% bc_R_fun = @(U_n,U_o,K,H,p)(U_n(1,end)-U_n(1,end-1));
% bc_L_jac_fun = @(U_n,U_o,K,H,p)[1,-1] ;
% bc_R_jac_fun = @(U_n,U_o,K,H,p)[-1,1];

% -------------------------------------------------------------------------
% Finite Difference Code
% -------------------------------------------------------------------------

p.something = 0;


U_n = sqrt(2)*sech(dom);
U_o = U_n-0.05*exp(-10*dom.^2);

% time steps
time_steps = 100;

K = 0.01;

% delta x
H = dom(2)-dom(1);

% time evolution code tolerance
tol = 1e-8;

% figure;
U_orig = U_o;
Ubar = U_n;
for j = 1:time_steps
    
    j
    
    U_n = finite_diff_advance(U_n,U_o,K,H,p,tol,@fd_F,@fd_jac, ...
    bc_L_fun,bc_R_fun,bc_L_jac_fun,bc_R_jac_fun);

    clf;
    hold on;
    plot(dom,Ubar,'--r','LineWidth',4);
    plot(dom,U_orig,'-k','LineWidth',1);
    plot(dom,U_n,'-g','LineWidth',2);
    drawnow;
    pause(0.1);

    U_o = U_n;

end
