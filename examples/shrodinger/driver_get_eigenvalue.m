beep off; clc; clear all; close all;
addpath(strcat(pwd,'/../../core/matlab'))

%
% parameters
%

p.a = 3; % a > 0
p.gamma = 7; % gamma >= 1

%
% dependent variables
%

p.mu = (p.a-sqrt(p.gamma^2-1))/(p.a+sqrt(p.gamma^2-1));
p.nu = 1/(p.a+sqrt(p.gamma^2-1));

%
% structure variables
%

s.I=300; 
s.R=s.I;
s.L=-s.I;


% lambda = 0.8;
% M = A(s.R,lambda,s,p);
% eig(M)
% return


% lambda = 10;
% 
% Q = 0.5*sqrt((p.mu+1)^2-4*(lambda^2+2*lambda*p.nu+p.mu));
% 
% C = 0.5*(1+p.mu);
% 
% mu1 = (sqrt(Q+C));
% mu2 = (sqrt(C-Q));
% 
% PL = [mu1,-1,0,0;
%       0, 0, mu2, -1];
%   
% PR = [-mu1,-1,0,0;
%       0, 0, -mu2,-1];
% 
% eig(A(s.L,lambda,s,p))

% return

% solve for the eigenfunction

% solve for the eigenfunction
%

L = s.I;
eigfun_ode = @(x,y,lambda)([A(x,lambda,s,p)*y(1:4); ...
    A(x+L,lambda,s,p)*y(5:8)]);



fun_bc = @(ya,yb,lambda)(bc_eigfun(ya,yb,lambda,p,s));

guess_fun = @(x)[sech(x);-sech(x)*tanh(x);sech(x);-sech(x)*tanh(x)];

lambda0 = 0.0531 + 1.9312i;

guess = @(x)[guess_fun(x);guess_fun(x+L)];

solinit = bvpinit(linspace(-L,0,30),guess,lambda0);
options = bvpset('RelTol',1e-10,'AbsTol',1e-10,'Nmax', 20000);
eigfun_sol = bvp5c(eigfun_ode,fun_bc,solinit,options);
lam = eigfun_sol.parameters;

%% plot eigenfunction
x = linspace(-L,0,500);
y = deval(eigfun_sol,x);
hold on;
plot(x,real(y(1:4,:)),'-g','LineWidth',2)
plot(x+L,real(y(5:8,:)),'-b','LineWidth',2)

return


return

%
% structure variables
%

[s,e,m,c] = emcset(s,'front',[2,2],'default'); % default for capillarity is reg_reg_polar






% m.method=@drury_no_radial;
% refine the contour to achieve desired relative tolerance
c.refine = 'on';
% Set the desired relative tolerance of the Evans function output
c.tol = 0.6;
% Set the number of Kato steps between the initial points
c.ksteps = 2^8;
c.lambda_steps = 0;
c.pic_stats = 'on';
c.root_fun = @contour;
% 
tol = 0.01;%10^(-2);
% 
beg_time = tic;
box = [0,1.8, 0.1, 2];
rt1 = root_solver1(box,tol,p,s,e,m,c);
% 
% box = [0.7,-2.5, 1];
% rt2 = root_solver2(box,tol,p,s,e,m,c);
% % 
return



% return

% display a waitbar
c.stats = 'print';
c.refine = 'on';
 
% m.ode_fun = @ode15s;

%
% preimage contour
%




% c.lambda_steps = 0;
% c.ksteps = 2^5;
% pnts = 30;
% preimage = linspace(3,0.7,pnts+(pnts-1)*c.ksteps);
% %
% % compute Evans function
% %
% tic
% [halfw,dom]=contour(c,s,p,m,e,preimage);
% time = toc
% plot(dom,halfw,'-k','LineWidth',2);
% return


% c.ksteps = 2^10;
circpnts=30; imagpnts=30; innerpnts = 10; r=2; spread=2; zerodist = 1.8;
preimage=semicirc2(circpnts,imagpnts,innerpnts,c.ksteps,r,spread,zerodist,c.lambda_steps)+0.1;
tic
[halfw,dom]=contour(c,s,p,m,e,preimage);
halfw = halfw/halfw(1);
time = toc
w = [halfw fliplr(conj(halfw))];

% 
% process and display data
%

wnd=winding_number(w);
fprintf('Winding Number: %1.1d\n',wnd);
% plot(w/w(end),'.-k');
plot_evans(w);


