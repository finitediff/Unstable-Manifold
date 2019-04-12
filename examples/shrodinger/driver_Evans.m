
addpath(strcat(pwd,'/../../core/matlab'))
beep off; clc; clear all; close all;

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

s.I=40; 
s.R=s.I;
s.L=-s.I;


% lambda = 0.8;
% M = A(s.R,lambda,s,p);
% eig(M)
% return


lambda = 10;

Q = 0.5*sqrt((p.mu+1)^2-4*(lambda^2+2*lambda*p.nu+p.mu));

C = 0.5*(1+p.mu);

mu1 = (sqrt(Q+C));
mu2 = (sqrt(C-Q));

PL = [mu1,-1,0,0;
      0, 0, mu2, -1];
  
PR = [-mu1,-1,0,0;
      0, 0, -mu2,-1];

eig(A(s.L,lambda,s,p))

% return

% solve for the eigenfunction

% solve for the eigenfunction
%%

L = s.I;
eigfun_ode = @(x,y,lambda)([A(x,lambda,s,p)*y(1:4); ...
    A(x+L,lambda,s,p)*y(5:8)]);



fun_bc = @(ya,yb,lambda)(bc_eigfun(ya,yb,lambda,p,s));

guess_fun = @(x)[sech(x);-sech(x)*tanh(x);sech(x);-sech(x)*tanh(x)];
lambda0 = 2+0.5*1i;

guess = @(x)[guess_fun(x);guess_fun(x+L)];

solinit = bvpinit(linspace(-L,0,30),guess,lambda0);
options = bvpset('RelTol',1e-6,'AbsTol',1e-6,'Nmax', 20000);
eigfun_sol = bvp5c(eigfun_ode,fun_bc,solinit,options);
lam = eigfun_sol.parameters;

%% plot eigenfunction
x = linspace(-L,0,500);
y = deval(eigfun_sol,x);
hold on;
plot(x,imag(y(1:4,:)),'-g','LineWidth',2)
plot(x+L,imag(y(5:8,:)),'-b','LineWidth',2)

return


return

%
% structure variables
%

[s,e,m,c] = emcset(s,'front',[2,2],'default'); % default for capillarity is reg_reg_polar



% display a waitbar
c.stats = 'print';
c.refine = 'on';
 
% m.ode_fun = @ode15s;

%
% preimage contour
%


circpnts=30; imagpnts=30; innerpnts = 10; r=5; spread=2; zerodist = 0.7;
preimage=semicirc2(circpnts,imagpnts,innerpnts,c.ksteps,r,spread,zerodist,c.lambda_steps);


%
% compute Evans function
%

tic
halfw=contour(c,s,p,m,e,preimage);
time = toc
w = [halfw fliplr(conj(halfw))];

% 
% process and display data
%

wnd=winding_number(w);
fprintf('Winding Number: %1.1d\n',wnd);
% plot(w/w(end),'.-k');
plot_evans(w);


