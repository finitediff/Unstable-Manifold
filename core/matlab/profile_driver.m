clear all; close all; clc; beep off; curr_dir = cd;

%
% parameters
%

p.gamma = 1/9;
p.alpha = 20;

% controls
s.I = 30; % spatial infinity

d.p = p;

% dependent controls
s.L = -s.I;
s.R = s.I;

% solve the profile
ode_fun = @(x,y)(profile_ode(x,y,p));
bc_fun = @(ya,yb)(profile_bc(ya,yb,p));
x = linspace(s.L,0,30);
guess = @(x)(profile_guess(x,p)); % NOT USING CONTINUATION
solinit = bvpinit(x,guess);
options = bvpset('RelTol',1e-6,'AbsTol',1e-6);
lastwarn('');
s.sol = bvp5c(ode_fun,bc_fun,solinit,options);

plot(s.sol.x,s.sol.y)

